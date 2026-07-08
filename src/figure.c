#include "field_image.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define PI 3.141592653589793238462643383279502884
#if defined(__GNUC__)
#define UNUSED_FUNCTION __attribute__((unused))
#else
#define UNUSED_FUNCTION
#endif

typedef struct {
    int width;
    int height;
    unsigned char *rgb;
} Canvas;

typedef struct {
    int width;
    int height;
    double *data;
    double min;
    double max;
} Field;

typedef struct {
    double *theta;
    double *confidence;
} Orientation;

typedef struct {
    double x;
    double y;
    double theta;
    double bias;
} SceneGrain;

typedef struct {
    int width;
    int height;
    int count;
    SceneGrain *grains;
} Scene;

static int compare_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t bytes) {
    void *ptr = malloc(bytes);
    if (ptr == NULL) die("Out of memory.");
    return ptr;
}

static int wrap_index(int value, int length) {
    value %= length;
    return value < 0 ? value + length : value;
}

static void canvas_init(Canvas *canvas, int width, int height) {
    canvas->width = width;
    canvas->height = height;
    canvas->rgb = xmalloc((size_t)width * height * 3);
    memset(canvas->rgb, 255, (size_t)width * height * 3);
}

static void canvas_free(Canvas *canvas) {
    free(canvas->rgb);
    memset(canvas, 0, sizeof(*canvas));
}

static void pixel(Canvas *canvas, int x, int y, int r, int g, int b) {
    if (x < 0 || y < 0 || x >= canvas->width || y >= canvas->height) return;
    unsigned char *p = canvas->rgb + ((size_t)y * canvas->width + x) * 3;
    p[0] = (unsigned char)r;
    p[1] = (unsigned char)g;
    p[2] = (unsigned char)b;
}

static void blend_pixel(Canvas *canvas, int x, int y, int r, int g, int b,
                        double alpha) {
    if (x < 0 || y < 0 || x >= canvas->width || y >= canvas->height) return;
    if (alpha <= 0.0) return;
    if (alpha > 1.0) alpha = 1.0;
    unsigned char *p = canvas->rgb + ((size_t)y * canvas->width + x) * 3;
    p[0] = (unsigned char)(p[0] * (1.0 - alpha) + r * alpha + 0.5);
    p[1] = (unsigned char)(p[1] * (1.0 - alpha) + g * alpha + 0.5);
    p[2] = (unsigned char)(p[2] * (1.0 - alpha) + b * alpha + 0.5);
}

static void soft_point(Canvas *canvas, double x, double y, int r, int g, int b,
                       double alpha) {
    int ix = (int)floor(x);
    int iy = (int)floor(y);
    for (int oy = -1; oy <= 1; ++oy) {
        for (int ox = -1; ox <= 1; ++ox) {
            double dx = x - (ix + ox);
            double dy = y - (iy + oy);
            double w = exp(-(dx * dx + dy * dy) / 0.42);
            blend_pixel(canvas, ix + ox, iy + oy, r, g, b, alpha * w);
        }
    }
}

static void UNUSED_FUNCTION soft_line(Canvas *canvas, double x0, double y0,
                                      double x1, double y1, int r, int g,
                                      int b, double alpha) {
    double dx = x1 - x0;
    double dy = y1 - y0;
    int steps = (int)(sqrt(dx * dx + dy * dy) * 1.8) + 1;
    for (int i = 0; i <= steps; ++i) {
        double t = steps == 0 ? 0.0 : (double)i / steps;
        soft_point(canvas, x0 + dx * t, y0 + dy * t, r, g, b, alpha);
    }
}

static void UNUSED_FUNCTION soft_line_rect(Canvas *canvas, double x0, double y0,
                                           double x1, double y1, int left,
                                           int top, int width, int height,
                                           int r, int g, int b,
                                           double alpha) {
    double dx = x1 - x0;
    double dy = y1 - y0;
    int steps = (int)(sqrt(dx * dx + dy * dy) * 1.8) + 1;
    for (int i = 0; i <= steps; ++i) {
        double t = steps == 0 ? 0.0 : (double)i / steps;
        double x = x0 + dx * t;
        double y = y0 + dy * t;
        if (x >= left && x <= left + width && y >= top && y <= top + height) {
            soft_point(canvas, x, y, r, g, b, alpha);
        }
    }
}

static void rect(Canvas *canvas, int x, int y, int w, int h, int r, int g,
                 int b) {
    for (int yy = y; yy < y + h; ++yy) {
        for (int xx = x; xx < x + w; ++xx) pixel(canvas, xx, yy, r, g, b);
    }
}

static void disc(Canvas *canvas, int cx, int cy, int radius, int r, int g,
                 int b) {
    for (int y = cy - radius; y <= cy + radius; ++y) {
        for (int x = cx - radius; x <= cx + radius; ++x) {
            int dx = x - cx;
            int dy = y - cy;
            if (dx * dx + dy * dy <= radius * radius) {
                pixel(canvas, x, y, r, g, b);
            }
        }
    }
}

static void line(Canvas *canvas, int x0, int y0, int x1, int y1, int r, int g,
                 int b) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        pixel(canvas, x0, y0, r, g, b);
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
}

static void triangle(Canvas *canvas, int x0, int y0, int x1, int y1, int x2,
                     int y2, int r, int g, int b) {
    int min_x = x0 < x1 ? x0 : x1;
    if (x2 < min_x) min_x = x2;
    int max_x = x0 > x1 ? x0 : x1;
    if (x2 > max_x) max_x = x2;
    int min_y = y0 < y1 ? y0 : y1;
    if (y2 < min_y) min_y = y2;
    int max_y = y0 > y1 ? y0 : y1;
    if (y2 > max_y) max_y = y2;
    if (min_x < 0) min_x = 0;
    if (min_y < 0) min_y = 0;
    if (max_x >= canvas->width) max_x = canvas->width - 1;
    if (max_y >= canvas->height) max_y = canvas->height - 1;
    int area = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
    if (area == 0) return;
    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            int w0 = (x1 - x) * (y2 - y) - (y1 - y) * (x2 - x);
            int w1 = (x2 - x) * (y0 - y) - (y2 - y) * (x0 - x);
            int w2 = (x0 - x) * (y1 - y) - (y0 - y) * (x1 - x);
            if ((area > 0 && w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                (area < 0 && w0 <= 0 && w1 <= 0 && w2 <= 0)) {
                pixel(canvas, x, y, r, g, b);
            }
        }
    }
}

#ifndef _WIN32
static const unsigned char *glyph(char c) {
    static const unsigned char blank[7] = {0, 0, 0, 0, 0, 0, 0};
    static const unsigned char qmark[7] = {14, 17, 1, 2, 4, 0, 4};
    static const unsigned char glyphs[43][7] = {
        {14, 17, 19, 21, 25, 17, 14}, {4, 12, 4, 4, 4, 4, 14},
        {14, 17, 1, 2, 4, 8, 31},    {30, 1, 1, 14, 1, 1, 30},
        {2, 6, 10, 18, 31, 2, 2},    {31, 16, 16, 30, 1, 1, 30},
        {14, 16, 16, 30, 17, 17, 14}, {31, 1, 2, 4, 8, 8, 8},
        {14, 17, 17, 14, 17, 17, 14}, {14, 17, 17, 15, 1, 1, 14},
        {14, 17, 17, 31, 17, 17, 17}, {30, 17, 17, 30, 17, 17, 30},
        {14, 17, 16, 16, 16, 17, 14}, {30, 17, 17, 17, 17, 17, 30},
        {31, 16, 16, 30, 16, 16, 31}, {31, 16, 16, 30, 16, 16, 16},
        {14, 17, 16, 23, 17, 17, 15}, {17, 17, 17, 31, 17, 17, 17},
        {14, 4, 4, 4, 4, 4, 14},      {7, 2, 2, 2, 18, 18, 12},
        {17, 18, 20, 24, 20, 18, 17}, {16, 16, 16, 16, 16, 16, 31},
        {17, 27, 21, 21, 17, 17, 17}, {17, 25, 21, 19, 17, 17, 17},
        {14, 17, 17, 17, 17, 17, 14}, {30, 17, 17, 30, 16, 16, 16},
        {14, 17, 17, 17, 21, 18, 13}, {30, 17, 17, 30, 20, 18, 17},
        {15, 16, 16, 14, 1, 1, 30},   {31, 4, 4, 4, 4, 4, 4},
        {17, 17, 17, 17, 17, 17, 14}, {17, 17, 17, 17, 10, 10, 4},
        {17, 17, 17, 21, 21, 21, 10}, {17, 17, 10, 4, 10, 17, 17},
        {17, 17, 10, 4, 4, 4, 4},     {31, 1, 2, 4, 8, 16, 31},
        {0, 0, 0, 0, 0, 0, 0},        {0, 4, 4, 31, 4, 4, 0},
        {0, 0, 0, 31, 0, 0, 0},       {0, 0, 0, 0, 0, 12, 12},
        {2, 2, 4, 4, 8, 8, 16},       {14, 16, 16, 16, 16, 16, 14},
        {14, 1, 1, 1, 1, 1, 14}};
    if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    if (c >= '0' && c <= '9') return glyphs[c - '0'];
    if (c >= 'A' && c <= 'Z') return glyphs[10 + c - 'A'];
    if (c == ' ') return glyphs[36];
    if (c == '+') return glyphs[37];
    if (c == '-') return glyphs[38];
    if (c == '.') return glyphs[39];
    if (c == '/') return glyphs[40];
    if (c == '(') return glyphs[41];
    if (c == ')') return glyphs[42];
    if (c == '\'' || c == '_') return glyphs[38];
    return c == '\0' ? blank : qmark;
}
#endif

#ifdef _WIN32
static void text_w(Canvas *canvas, int x, int y, const wchar_t *s, int scale,
                   int r, int g, int b) {
    if (s == NULL || *s == L'\0') return;

    BITMAPINFO info;
    memset(&info, 0, sizeof(info));
    info.bmiHeader.biSize = sizeof(info.bmiHeader);
    info.bmiHeader.biWidth = canvas->width;
    info.bmiHeader.biHeight = -canvas->height;
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;

    void *bits = NULL;
    HDC screen = GetDC(NULL);
    HDC dc = CreateCompatibleDC(screen);
    HBITMAP bitmap = CreateDIBSection(screen, &info, DIB_RGB_COLORS, &bits,
                                      NULL, 0);
    ReleaseDC(NULL, screen);
    if (dc == NULL || bitmap == NULL || bits == NULL) {
        if (bitmap != NULL) DeleteObject(bitmap);
        if (dc != NULL) DeleteDC(dc);
        return;
    }

    HGDIOBJ old_bitmap = SelectObject(dc, bitmap);
    unsigned char *dib = bits;
    for (int yy = 0; yy < canvas->height; ++yy) {
        for (int xx = 0; xx < canvas->width; ++xx) {
            size_t src = ((size_t)yy * canvas->width + xx) * 3;
            size_t dst = ((size_t)yy * canvas->width + xx) * 4;
            dib[dst + 0] = canvas->rgb[src + 2];
            dib[dst + 1] = canvas->rgb[src + 1];
            dib[dst + 2] = canvas->rgb[src + 0];
            dib[dst + 3] = 0;
        }
    }

    int font_height = scale <= 1 ? 17 : scale * 13;
    int weight = scale >= 4 ? FW_NORMAL : FW_REGULAR;
    HFONT font = CreateFontW(-font_height, 0, 0, 0, weight, FALSE, FALSE,
                             FALSE, DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
                             CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                             DEFAULT_PITCH | FF_SWISS, L"Arial");
    HGDIOBJ old_font = SelectObject(dc, font);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, RGB(r, g, b));
    TextOutW(dc, x, y, s, (int)wcslen(s));
    GdiFlush();

    for (int yy = 0; yy < canvas->height; ++yy) {
        for (int xx = 0; xx < canvas->width; ++xx) {
            size_t src = ((size_t)yy * canvas->width + xx) * 4;
            size_t dst = ((size_t)yy * canvas->width + xx) * 3;
            canvas->rgb[dst + 0] = dib[src + 2];
            canvas->rgb[dst + 1] = dib[src + 1];
            canvas->rgb[dst + 2] = dib[src + 0];
        }
    }

    SelectObject(dc, old_font);
    SelectObject(dc, old_bitmap);
    DeleteObject(font);
    DeleteObject(bitmap);
    DeleteDC(dc);
}

static void text(Canvas *canvas, int x, int y, const char *s, int scale, int r,
                 int g, int b) {
    wchar_t wide[256];
    int i = 0;
    for (; s[i] != '\0' && i < (int)(sizeof(wide) / sizeof(wide[0])) - 1; ++i) {
        wide[i] = (unsigned char)s[i];
    }
    wide[i] = L'\0';
    text_w(canvas, x, y, wide, scale, r, g, b);
}
#else
static void text(Canvas *canvas, int x, int y, const char *s, int scale, int r,
                 int g, int b) {
    int cursor = x;
    for (; *s; ++s) {
        const unsigned char *rows = glyph(*s);
        for (int yy = 0; yy < 7; ++yy) {
            for (int xx = 0; xx < 5; ++xx) {
                if (rows[yy] & (1 << (4 - xx))) {
                    rect(canvas, cursor + xx * scale, y + yy * scale, scale,
                         scale, r, g, b);
                }
            }
        }
        cursor += 6 * scale;
    }
}
#endif

#ifndef _WIN32
static void symbol(Canvas *canvas, int x, int y, const unsigned char rows[7],
                   int scale, int r, int g, int b) {
    for (int yy = 0; yy < 7; ++yy) {
        for (int xx = 0; xx < 5; ++xx) {
            if (rows[yy] & (1 << (4 - xx))) {
                rect(canvas, x + xx * scale, y + yy * scale, scale, scale, r,
                     g, b);
            }
        }
    }
}

static void symbol_pi(Canvas *canvas, int x, int y, int scale, int r, int g,
                      int b) {
    static const unsigned char rows[7] = {31, 17, 17, 17, 17, 17, 17};
    symbol(canvas, x, y, rows, scale, r, g, b);
}

static void symbol_sigma(Canvas *canvas, int x, int y, int scale, int r, int g,
                         int b) {
    static const unsigned char rows[7] = {0, 0, 15, 18, 18, 18, 12};
    symbol(canvas, x, y, rows, scale, r, g, b);
}
#endif

static void angle_label(Canvas *canvas, int x, int y, int negative, int scale,
                        int r, int g, int b) {
#ifdef _WIN32
    text_w(canvas, x, y, negative ? L"-\x03c0/12" : L"\x03c0/12", scale, r, g,
           b);
#else
    int cursor = x;
    if (negative) {
        text(canvas, cursor, y, "-", scale, r, g, b);
        cursor += 6 * scale;
    }
    symbol_pi(canvas, cursor, y, scale, r, g, b);
    cursor += 6 * scale;
    text(canvas, cursor, y, "/12", scale, r, g, b);
#endif
}

static void sigma_x_label(Canvas *canvas, int x, int y, int scale, int r,
                          int g, int b) {
#ifdef _WIN32
    text_w(canvas, x, y, L"\x03c3", scale, r, g, b);
    text(canvas, x + 11 * scale, y + 7 * scale, "x", 1, r, g, b);
#else
    symbol_sigma(canvas, x, y, scale, r, g, b);
    text(canvas, x + 6 * scale, y + 4 * scale, "X", 1, r, g, b);
#endif
}

static Field read_field(const char *filename, int width, int height) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Missing field file: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    Field field = {width, height, xmalloc((size_t)width * height * sizeof(double)),
                   DBL_MAX, -DBL_MAX};
    for (int i = 0; i < width * height; ++i) {
        if (fscanf(file, "%lf", &field.data[i]) != 1) {
            fprintf(stderr, "Invalid field data in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        if (field.data[i] < field.min) field.min = field.data[i];
        if (field.data[i] > field.max) field.max = field.data[i];
    }
    fclose(file);
    if (field.max == field.min) field.max = field.min + 1.0;
    return field;
}

static void free_field(Field *field) {
    free(field->data);
    memset(field, 0, sizeof(*field));
}

static double sample_field(const Field *field, double x, double y) {
    if (x < 0.0) x = 0.0;
    if (y < 0.0) y = 0.0;
    if (x > field->width - 1.0) x = field->width - 1.0;
    if (y > field->height - 1.0) y = field->height - 1.0;
    int x0 = (int)x, y0 = (int)y;
    int x1 = x0 + 1 < field->width ? x0 + 1 : x0;
    int y1 = y0 + 1 < field->height ? y0 + 1 : y0;
    double fx = x - x0, fy = y - y0;
    double v00 = field->data[(size_t)y0 * field->width + x0];
    double v10 = field->data[(size_t)y0 * field->width + x1];
    double v01 = field->data[(size_t)y1 * field->width + x0];
    double v11 = field->data[(size_t)y1 * field->width + x1];
    return (v00 + fx * (v10 - v00)) * (1.0 - fy) +
           (v01 + fx * (v11 - v01)) * fy;
}

static void hsv(double h, double s, double v, int *r, int *g, int *b) {
    h = h - floor(h);
    double c = v * s;
    double hp = h * 6.0;
    double x = c * (1.0 - fabs(fmod(hp, 2.0) - 1.0));
    double rr = 0, gg = 0, bb = 0;
    if (hp < 1) {
        rr = c; gg = x;
    } else if (hp < 2) {
        rr = x; gg = c;
    } else if (hp < 3) {
        gg = c; bb = x;
    } else if (hp < 4) {
        gg = x; bb = c;
    } else if (hp < 5) {
        rr = x; bb = c;
    } else {
        rr = c; bb = x;
    }
    double m = v - c;
    *r = (int)((rr + m) * 255.0 + 0.5);
    *g = (int)((gg + m) * 255.0 + 0.5);
    *b = (int)((bb + m) * 255.0 + 0.5);
}

static double wrap_theta(double theta) {
    double period = PI / 3.0;
    theta = fmod(theta + PI / 6.0, period);
    if (theta < 0.0) theta += period;
    return theta - PI / 6.0;
}

static void orientation_rgb(double theta, double confidence, int *r, int *g,
                            int *b);
static void render_colorbar(Canvas *canvas, int x, int y, int w, int h);
static void sigma_x_label(Canvas *canvas, int x, int y, int scale, int r,
                          int g, int b);

static double scene_rand(unsigned int *state) {
    *state = *state * 1664525u + 1013904223u;
    return (double)(*state >> 8) / 16777215.0;
}

static Scene UNUSED_FUNCTION scene_create(int width, int height, int count) {
    Scene scene = {width, height, count,
                   xmalloc((size_t)count * sizeof(*scene.grains))};
    int cols = (int)ceil(sqrt((double)count * width / height));
    int rows = (count + cols - 1) / cols;
    unsigned int state = 20260707u;
    int index = 0;
    for (int row = 0; row < rows && index < count; ++row) {
        for (int col = 0; col < cols && index < count; ++col) {
            double cell_w = (double)width / cols;
            double cell_h = (double)height / rows;
            scene.grains[index].x =
                (col + 0.5 + 0.70 * (scene_rand(&state) - 0.5)) * cell_w;
            scene.grains[index].y =
                (row + 0.5 + 0.70 * (scene_rand(&state) - 0.5)) * cell_h;
            scene.grains[index].theta =
                -PI / 12.0 + scene_rand(&state) * PI / 6.0;
            scene.grains[index].bias = 0.82 + 0.34 * scene_rand(&state);
            ++index;
        }
    }
    return scene;
}

static void UNUSED_FUNCTION scene_free(Scene *scene) {
    free(scene->grains);
    memset(scene, 0, sizeof(*scene));
}

static double scene_delta(double a, double b, double length) {
    double d = a - b;
    if (d > 0.5 * length) d -= length;
    if (d < -0.5 * length) d += length;
    return d;
}

static int scene_nearest(const Scene *scene, double x, double y, double *dist,
                         double *second_dist) {
    int best = 0;
    double best_d2 = HUGE_VAL;
    double second_d2 = HUGE_VAL;
    double qx = x + 1.4 * sin(0.033 * y) + 0.6 * sin(0.027 * x + 1.7);
    double qy = y + 0.9 * sin(0.039 * x - 0.4) + 0.5 * sin(0.029 * y);
    for (int i = 0; i < scene->count; ++i) {
        double dx = scene_delta(qx, scene->grains[i].x, scene->width);
        double dy = scene_delta(qy, scene->grains[i].y, scene->height);
        double d2 = dx * dx + 2.85 * dy * dy;
        if (d2 < best_d2) {
            second_d2 = best_d2;
            best_d2 = d2;
            best = i;
        } else if (d2 < second_d2) {
            second_d2 = d2;
        }
    }
    if (dist != NULL) *dist = sqrt(best_d2);
    if (second_dist != NULL) *second_dist = sqrt(second_d2);
    return best;
}

static int scene_boundary(const Scene *scene, double x, double y, int id) {
    int right = scene_nearest(scene, x + 0.9, y, NULL, NULL);
    int down = scene_nearest(scene, x, y + 0.9, NULL, NULL);
    return right != id || down != id;
}

static void scene_color(const Scene *scene, int id, int *r, int *g, int *b) {
    orientation_rgb(scene->grains[id].theta, 0.95, r, g, b);
}

static int scene_gray_texture(const Scene *scene, int id, double x, double y,
                              double edge) {
    double theta = scene->grains[id].theta;
    double c = cos(theta);
    double s = sin(theta);
    double rx = scene_delta(x, scene->grains[id].x, scene->width);
    double ry = scene_delta(y, scene->grains[id].y, scene->height);
    double u = c * rx + s * ry;
    double v = -s * rx + c * ry;
    double l1 = fabs(sin(1.62 * u));
    double l2 = fabs(sin(1.62 * (0.50 * u + 0.866 * v)));
    double l3 = fabs(sin(1.62 * (-0.50 * u + 0.866 * v)));
    double line_strength = fmax(fmax(1.0 - l1, 1.0 - l2), 1.0 - l3);
    if (line_strength < 0.0) line_strength = 0.0;
    if (line_strength > 1.0) line_strength = 1.0;
    double dots = 0.5 + 0.5 * sin(1.62 * u) * sin(1.62 * v);
    int gray = 188 - (int)(58.0 * line_strength + 18.0 * dots);
    if (edge < 1.8) gray = 176 + (int)(8.0 * sin(0.55 * x + 0.37 * y));
    if (gray < 94) gray = 94;
    if (gray > 218) gray = 218;
    return gray;
}

static double scene_growth_radius(const Scene *scene, int id, double x,
                                  double y, double base_radius, int frame) {
    double dx = scene_delta(x, scene->grains[id].x, scene->width);
    double dy = scene_delta(y, scene->grains[id].y, scene->height);
    double angle = atan2(dy, dx);
    double amount = frame == 1 ? 0.30 : (frame == 2 ? 0.55 : 0.42);
    double wobble = 1.0 + amount * (0.045 * sin(3.0 * angle + 0.3 * id) +
                                    0.030 * sin(5.0 * angle + 0.13 * id) +
                                    0.018 * sin(0.045 * x - 0.038 * y + id));
    return base_radius * scene->grains[id].bias * wobble;
}

static void UNUSED_FUNCTION render_scene_growth(Canvas *canvas,
                                                const Scene *scene, int x0,
                                                int y0, int w, int h,
                                                int frame) {
    static const double radii[] = {3.2, 11.5, 18.0, 38.0};
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    for (int y = 0; y < h; ++y) {
        double sy = (double)y / (h - 1) * (scene->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (scene->width - 1);
            double dist, second;
            int id = scene_nearest(scene, sx, sy, &dist, &second);
            double radius =
                scene_growth_radius(scene, id, sx, sy, radii[frame], frame);
            double edge = radius - dist;
            if (frame < 3 && edge < -2.5) continue;
            int gray;
            if (frame == 0) {
                double spot = exp(-(dist * dist) / (2.0 * radius * radius));
                if (spot < 0.12) continue;
                gray = 248 - (int)(92.0 * spot);
            } else {
                gray = scene_gray_texture(scene, id, sx, sy, edge);
                if (edge < 0.0) {
                    double t = (edge + 2.5) / 2.5;
                    int halo = 250 - (int)(22.0 * t);
                    gray = (int)(halo * (1.0 - t) + gray * t);
                }
                if (frame == 3 &&
                    fabs(second - dist + 0.7 * sin(0.08 * sx + 0.05 * sy)) <
                        0.48) {
                    gray = 188;
                }
            }
            pixel(canvas, x0 + x, y0 + y, gray, gray, gray);
        }
    }
}

static void UNUSED_FUNCTION render_scene_orientation(Canvas *canvas,
                                                     const Scene *scene, int x0,
                                                     int y0, int w, int h) {
    for (int y = 0; y < h; ++y) {
        double sy = (double)y / (h - 1) * (scene->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (scene->width - 1);
            int id = scene_nearest(scene, sx, sy, NULL, NULL);
            int r, g, b;
            scene_color(scene, id, &r, &g, &b);
            if (scene_boundary(scene, sx, sy, id)) {
                r = (r * 9 + 235) / 10;
                g = (g * 9 + 235) / 10;
                b = (b * 9 + 235) / 10;
            }
            pixel(canvas, x0 + x, y0 + y, r, g, b);
        }
    }
}

static void UNUSED_FUNCTION render_scene_3d(Canvas *canvas, const Scene *scene,
                                            int x0, int y0) {
    render_colorbar(canvas, x0 + 142, y0 - 34, 120, 16);
    int grid_w = 47;
    int grid_h = 31;
    int *px = xmalloc((size_t)grid_w * grid_h * sizeof(*px));
    int *py = xmalloc((size_t)grid_w * grid_h * sizeof(*py));
    int *rr = xmalloc((size_t)grid_w * grid_h * sizeof(*rr));
    int *gg = xmalloc((size_t)grid_w * grid_h * sizeof(*gg));
    int *bb = xmalloc((size_t)grid_w * grid_h * sizeof(*bb));
    for (int gy = 0; gy < grid_h; ++gy) {
        double sy = (double)gy / (grid_h - 1) * (scene->height - 1);
        for (int gx = 0; gx < grid_w; ++gx) {
            double sx = (double)gx / (grid_w - 1) * (scene->width - 1);
            int id = scene_nearest(scene, sx, sy, NULL, NULL);
            double u = sx - scene->width / 2.0;
            double v = sy - scene->height / 2.0;
            double z = 6.7 * sin(0.024 * sx + 0.6) *
                           cos(0.032 * sy - 0.3) +
                       2.1 * sin(0.034 * sx - 0.020 * sy);
            size_t idx = (size_t)gy * grid_w + gx;
            px[idx] = x0 + 205 + (int)(0.95 * u - 0.38 * v);
            py[idx] = y0 + 145 + (int)(0.42 * v - z);
            scene_color(scene, id, &rr[idx], &gg[idx], &bb[idx]);
        }
    }
    for (int gy = grid_h - 2; gy >= 0; --gy) {
        for (int gx = 0; gx < grid_w - 1; ++gx) {
            size_t a = (size_t)gy * grid_w + gx;
            size_t b = a + 1;
            size_t c = a + grid_w;
            size_t d = c + 1;
            int r = (rr[a] + rr[b] + rr[c] + rr[d]) / 4;
            int g = (gg[a] + gg[b] + gg[c] + gg[d]) / 4;
            int bl = (bb[a] + bb[b] + bb[c] + bb[d]) / 4;
            r = (r * 7 + 245) / 8;
            g = (g * 7 + 245) / 8;
            bl = (bl * 7 + 245) / 8;
            triangle(canvas, px[a], py[a], px[b], py[b], px[c], py[c], r, g,
                     bl);
            triangle(canvas, px[b], py[b], px[d], py[d], px[c], py[c], r, g,
                     bl);
            if ((gx + gy) % 5 == 0) {
                line(canvas, px[a], py[a], px[b], py[b], r / 2, g / 2,
                     bl / 2);
            }
        }
    }
    free(px);
    free(py);
    free(rr);
    free(gg);
    free(bb);

    line(canvas, x0 + 75, y0 + 178, x0 + 362, y0 + 178, 160, 95, 40);
    line(canvas, x0 + 362, y0 + 178, x0 + 332, y0 + 203, 160, 95, 40);
    text(canvas, x0 + 205, y0 + 184, "L", 2, 160, 95, 40);
    line(canvas, x0 + 386, y0 + 98, x0 + 386, y0 + 178, 160, 95, 40);
    text(canvas, x0 + 392, y0 + 132, "L", 2, 160, 95, 40);
    for (int i = 0; i < 5; ++i) {
        int yy = y0 + 68 + i * 19;
        line(canvas, x0 + 4, yy, x0 + 48, yy, 0, 160, 220);
        line(canvas, x0 + 4, yy, x0 + 14, yy - 4, 0, 160, 220);
        line(canvas, x0 + 4, yy, x0 + 14, yy + 4, 0, 160, 220);
        line(canvas, x0 + 404, yy, x0 + 448, yy, 0, 160, 220);
        line(canvas, x0 + 448, yy, x0 + 438, yy - 4, 0, 160, 220);
        line(canvas, x0 + 448, yy, x0 + 438, yy + 4, 0, 160, 220);
    }
    sigma_x_label(canvas, x0 + 2, y0 + 42, 2, 0, 160, 220);
    sigma_x_label(canvas, x0 + 420, y0 + 42, 2, 0, 160, 220);
    line(canvas, x0 + 28, y0 + 203, x0 + 66, y0 + 203, 255, 40, 30);
    line(canvas, x0 + 28, y0 + 203, x0 + 28, y0 + 166, 40, 120, 255);
    line(canvas, x0 + 28, y0 + 203, x0 + 47, y0 + 184, 70, 190, 70);
    text(canvas, x0 + 69, y0 + 196, "x", 1, 255, 40, 30);
    text(canvas, x0 + 21, y0 + 151, "z", 1, 40, 120, 255);
    text(canvas, x0 + 50, y0 + 177, "y", 1, 70, 190, 70);
}

static void UNUSED_FUNCTION render_scene_boundaries(Canvas *canvas,
                                                    const Scene *scene, int x0,
                                                    int y0, int w, int h) {
    rect(canvas, x0, y0, w, h, 252, 252, 252);
    for (int i = 0; i <= 4; ++i) {
        int xx = x0 + i * w / 4;
        int yy = y0 + i * h / 4;
        line(canvas, xx, y0, xx, y0 + h, 172, 172, 172);
        line(canvas, x0, yy, x0 + w, yy, 172, 172, 172);
    }
    for (int y = 1; y < h - 1; ++y) {
            double sy = (double)y / (h - 1) * (scene->height - 1);
        for (int x = 1; x < w - 1; ++x) {
            double sx = (double)x / (w - 1) * (scene->width - 1);
            double curve = 0.65 * sin(0.08 * sx + 0.05 * sy);
            int id2 = scene_nearest(scene, sx + curve, sy - 0.45 * curve,
                                    NULL, NULL);
            if (scene_boundary(scene, sx + curve, sy - 0.45 * curve, id2) &&
                ((x + 3 * y) % 7 != 0)) {
                pixel(canvas, x0 + x, y0 + y, 232, 82, 72);
            }
        }
    }
    text(canvas, x0 + w / 2 - 7, y0 + h + 11, "N", 2, 40, 40, 50);
    text(canvas, x0 + w / 2 + 11, y0 + h + 26, "x", 1, 40, 40, 50);
    text(canvas, x0 - 34, y0 + h / 2 - 10, "N", 2, 40, 40, 50);
    text(canvas, x0 - 17, y0 + h / 2 + 5, "y", 1, 40, 40, 50);
    line(canvas, x0 - 18, y0 + h + 20, x0 + 22, y0 + h + 20, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h + 20, x0 - 18, y0 + h - 20, 50, 50, 50);
    line(canvas, x0 + 22, y0 + h + 20, x0 + 14, y0 + h + 16, 50, 50, 50);
    line(canvas, x0 + 22, y0 + h + 20, x0 + 14, y0 + h + 24, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h - 20, x0 - 22, y0 + h - 12, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h - 20, x0 - 14, y0 + h - 12, 50, 50, 50);
    text(canvas, x0 + 26, y0 + h + 14, "x", 1, 40, 40, 50);
    text(canvas, x0 - 26, y0 + h - 34, "y", 1, 40, 40, 50);
}

static Orientation compute_orientation(const Field *field, int radius) {
    Orientation out = {xmalloc((size_t)field->width * field->height *
                               sizeof(double)),
                       xmalloc((size_t)field->width * field->height *
                               sizeof(double))};
    const int block = 4;
    const int candidates = 37;
    const double dx = 0.7;
    const double dy = 0.7;
    const double lattice = 7.2552;
    const double qx = 2.0 * PI / lattice;
    const double qy = qx / sqrt(3.0);
    const double base[3][2] = {{qx, qy}, {qx, -qy}, {0.0, 2.0 * qy}};
    int bw = (field->width + block - 1) / block;
    int bh = (field->height + block - 1) / block;
    double *btheta = xmalloc((size_t)bw * bh * sizeof(double));
    double *bconf = xmalloc((size_t)bw * bh * sizeof(double));

    for (int by = 0; by < bh; ++by) {
        int cy = by * block + block / 2;
        if (cy >= field->height) cy = field->height - 1;
        for (int bx = 0; bx < bw; ++bx) {
            int cx = bx * block + block / 2;
            if (cx >= field->width) cx = field->width - 1;
            double mean = 0.0;
            int count = 0;
            for (int oy = -radius; oy <= radius; ++oy) {
                int yy = wrap_index(cy + oy, field->height);
                for (int ox = -radius; ox <= radius; ++ox) {
                    int xx = wrap_index(cx + ox, field->width);
                    mean += field->data[(size_t)yy * field->width + xx];
                    ++count;
                }
            }
            mean /= count;

            double best = -1.0;
            double second = -1.0;
            double best_theta = 0.0;
            for (int c = 0; c < candidates; ++c) {
                double theta = -PI / 6.0 +
                               (PI / 3.0) * c / (double)(candidates - 1);
                double cost = cos(theta);
                double sint = sin(theta);
                double score = 0.0;
                for (int k = 0; k < 3; ++k) {
                    double kx = base[k][0] * cost + base[k][1] * sint;
                    double ky = -base[k][0] * sint + base[k][1] * cost;
                    double re = 0.0;
                    double im = 0.0;
                    for (int oy = -radius; oy <= radius; ++oy) {
                        int yy = wrap_index(cy + oy, field->height);
                        double py = oy * dy;
                        for (int ox = -radius; ox <= radius; ++ox) {
                            int xx = wrap_index(cx + ox, field->width);
                            double px = ox * dx;
                            double value =
                                field->data[(size_t)yy * field->width + xx] -
                                mean;
                            double phase = kx * px + ky * py;
                            re += value * cos(phase);
                            im += value * sin(phase);
                        }
                    }
                    score += re * re + im * im;
                }
                if (score > best) {
                    second = best;
                    best = score;
                    best_theta = theta;
                } else if (score > second) {
                    second = score;
                }
            }
            btheta[(size_t)by * bw + bx] = wrap_theta(best_theta);
            bconf[(size_t)by * bw + bx] =
                best <= 0.0 ? 0.0 : (best - second) / best;
        }
    }

    for (int y = 0; y < field->height; ++y) {
        int by = y / block;
        for (int x = 0; x < field->width; ++x) {
            int bx = x / block;
            size_t dst = (size_t)y * field->width + x;
            size_t src = (size_t)by * bw + bx;
            out.theta[dst] = btheta[src];
            out.confidence[dst] = 0.55 + 0.45 * bconf[src];
        }
    }

    free(btheta);
    free(bconf);
    return out;
}

static void smooth_orientation(Orientation *orientation, int width, int height,
                               int iterations) {
    double *next_theta = xmalloc((size_t)width * height * sizeof(*next_theta));
    double *next_conf =
        xmalloc((size_t)width * height * sizeof(*next_conf));
    for (int iter = 0; iter < iterations; ++iter) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                size_t idx = (size_t)y * width + x;
                double center = orientation->theta[idx];
                double sum_c = 0.0;
                double sum_s = 0.0;
                double sum_conf = 0.0;
                int count = 0;
                for (int oy = -2; oy <= 2; ++oy) {
                    int yy = wrap_index(y + oy, height);
                    for (int ox = -2; ox <= 2; ++ox) {
                        int xx = wrap_index(x + ox, width);
                        size_t n = (size_t)yy * width + xx;
                        double theta = orientation->theta[n];
                        if (fabs(wrap_theta(theta - center)) > 0.070) continue;
                        double phase = 6.0 * theta;
                        double weight = orientation->confidence[n];
                        sum_c += weight * cos(phase);
                        sum_s += weight * sin(phase);
                        sum_conf += orientation->confidence[n];
                        ++count;
                    }
                }
                if (count == 0) {
                    next_theta[idx] = center;
                    next_conf[idx] = orientation->confidence[idx];
                } else {
                    next_theta[idx] = wrap_theta(atan2(sum_s, sum_c) / 6.0);
                    next_conf[idx] = sum_conf / count;
                }
            }
        }
        memcpy(orientation->theta, next_theta,
               (size_t)width * height * sizeof(*next_theta));
        memcpy(orientation->confidence, next_conf,
               (size_t)width * height * sizeof(*next_conf));
    }
    free(next_theta);
    free(next_conf);
}

static void free_orientation(Orientation *orientation) {
    free(orientation->theta);
    free(orientation->confidence);
    memset(orientation, 0, sizeof(*orientation));
}

static void orientation_rgb(double theta, double confidence, int *r, int *g,
                            int *b) {
    double limit = PI / 12.0;
    if (theta < -limit) theta = -limit;
    if (theta > limit) theta = limit;
    double t = (theta + limit) / (2.0 * limit);
    double hue = (2.0 / 3.0) * (1.0 - t);
    double sat = 0.70 + 0.25 * confidence;
    if (sat > 1.0) sat = 1.0;
    hsv(hue, sat, 0.90, r, g, b);
}

static void UNUSED_FUNCTION gray_range(const Field *fields, int count,
                                       double *low, double *high) {
    size_t total = 0;
    for (int i = 0; i < count; ++i) {
        total += (size_t)fields[i].width * fields[i].height;
    }
    double *sorted = xmalloc(total * sizeof(*sorted));
    size_t offset = 0;
    double minimum = HUGE_VAL;
    double maximum = -HUGE_VAL;
    for (int i = 0; i < count; ++i) {
        size_t n = (size_t)fields[i].width * fields[i].height;
        memcpy(sorted + offset, fields[i].data, n * sizeof(*sorted));
        offset += n;
        if (fields[i].min < minimum) minimum = fields[i].min;
        if (fields[i].max > maximum) maximum = fields[i].max;
    }
    qsort(sorted, total, sizeof(*sorted), compare_double);
    *low = sorted[total / 50];
    *high = sorted[total - total / 50 - 1];
    free(sorted);
    if (*high <= *low) {
        *low = minimum;
        *high = maximum;
    }
}

static void UNUSED_FUNCTION render_gray(Canvas *canvas, const Field *field,
                                        int x0, int y0, int w, int h,
                                        double low, double high, int frame) {
    double range = high - low;
    for (int y = 0; y < h; ++y) {
        double sy = field->height - 1.0 - (double)y / (h - 1) *
                                           (field->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (field->width - 1);
            double v = (sample_field(field, sx, sy) - low) / range;
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            double stripe = 0.012 * sin(0.42 * sx + 0.74 * sy) +
                            0.008 * sin(0.83 * sx - 0.34 * sy + 1.7);
            if (frame == 0) {
                double spot = 1.0 - v;
                if (spot < 0.43) {
                    v = 0.985;
                } else {
                    double t = (spot - 0.43) / 0.57;
                    if (t > 1.0) t = 1.0;
                    v = 0.97 - 0.44 * t;
                }
            } else if (frame < 3) {
                v = 0.94 - 0.62 * pow(1.0 - v + stripe, 0.82);
            } else {
                v = 0.88 - 0.48 * pow(1.0 - v + 0.35 * stripe, 0.90);
            }
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            int c = (int)(255.0 * v + 0.5);
            pixel(canvas, x0 + x, y0 + y, c, c, c);
        }
    }
}

static void UNUSED_FUNCTION render_orientation(Canvas *canvas,
                                               const Field *field,
                                               const Orientation *orientation,
                                               int x0, int y0, int w, int h,
                                               int draw_boundaries) {
    (void)field;
    for (int y = 0; y < h; ++y) {
        int sy = field->height - 1 -
                 (int)((double)y / (h - 1) * (field->height - 1) + 0.5);
        for (int x = 0; x < w; ++x) {
            double fx = (double)x / (w - 1) * (field->width - 1);
            int sx = (int)(fx + 0.5);
            int x0s = (int)floor(fx);
            int x1s = wrap_index(x0s + 1, field->width);
            double tx = fx - x0s;
            size_t idx = (size_t)sy * field->width + sx;
            size_t idx1 = (size_t)sy * field->width + x1s;
            double theta =
                atan2((1.0 - tx) * sin(6.0 * orientation->theta[idx]) +
                          tx * sin(6.0 * orientation->theta[idx1]),
                      (1.0 - tx) * cos(6.0 * orientation->theta[idx]) +
                          tx * cos(6.0 * orientation->theta[idx1])) /
                6.0;
            double confidence = (1.0 - tx) * orientation->confidence[idx] +
                                tx * orientation->confidence[idx1];
            int r, g, b;
            orientation_rgb(theta, confidence * 0.82, &r, &g, &b);
            r = (r * 3 + 245 * 2) / 5;
            g = (g * 3 + 245 * 2) / 5;
            b = (b * 3 + 245 * 2) / 5;
            if (draw_boundaries) {
                int sx1 = wrap_index(sx + 1, field->width);
                int sy1 = wrap_index(sy + 1, field->height);
                double dx = fabs(wrap_theta(orientation->theta[idx] -
                                            orientation->theta[(size_t)sy *
                                                field->width + sx1]));
                double dy = fabs(wrap_theta(orientation->theta[idx] -
                                            orientation->theta[(size_t)sy1 *
                                                field->width + sx]));
                if (dx > 0.18 || dy > 0.18) {
                    r = (r + 210) / 2;
                    g /= 3;
                    b /= 3;
                }
            }
            pixel(canvas, x0 + x, y0 + y, r, g, b);
        }
    }
}

static void render_colorbar(Canvas *canvas, int x, int y, int w, int h) {
    for (int xx = 0; xx < w; ++xx) {
        double t = (double)xx / (w - 1);
        int r, g, b;
        hsv((2.0 / 3.0) * (1.0 - t), 0.8, 0.95, &r, &g, &b);
        for (int yy = 0; yy < h; ++yy) pixel(canvas, x + xx, y + yy, r, g, b);
    }
    angle_label(canvas, x - 58, y + 1, 1, 1, 40, 40, 50);
    angle_label(canvas, x + w + 8, y + 1, 0, 1, 40, 40, 50);
}

static void render_vertical_colorbar(Canvas *canvas, int x, int y, int w,
                                     int h) {
    for (int yy = 0; yy < h; ++yy) {
        double t = 1.0 - (double)yy / (h - 1);
        int r, g, b;
        hsv((2.0 / 3.0) * (1.0 - t), 0.8, 0.95, &r, &g, &b);
        for (int xx = 0; xx < w; ++xx) pixel(canvas, x + xx, y + yy, r, g, b);
    }
    angle_label(canvas, x - 8, y - 20, 0, 1, 40, 40, 50);
    angle_label(canvas, x - 20, y + h + 8, 1, 1, 40, 40, 50);
}

static void render_plot(Canvas *canvas, int x, int y, int w, int h,
                        int variant) {
    (void)variant;
    rect(canvas, x, y, w, h, 250, 250, 250);
    line(canvas, x, y + h, x + w, y + h, 30, 30, 30);
    line(canvas, x, y, x, y + h, 30, 30, 30);
    line(canvas, x, y + h, x + w, y, 20, 20, 20);
    for (int value = 0; value <= 14; value += 2) {
        int xx = x + (int)(value / 14.0 * w + 0.5);
        int yy = y + h - (int)(value / 14.0 * h + 0.5);
        char label[8];
        snprintf(label, sizeof(label), "%d", value);
        line(canvas, xx, y + h, xx, y + h + 5, 30, 30, 30);
        line(canvas, x - 5, yy, x, yy, 30, 30, 30);
        text(canvas, xx - (value >= 10 ? 7 : 3), y + h + 8, label, 1, 20, 20,
             20);
        text(canvas, x - (value >= 10 ? 20 : 14), yy - 4, label, 1, 20, 20,
             20);
    }
    double dmin = 2.8 / 14.0;
    int dmin_y = y + h - (int)(dmin * h);
    line(canvas, x, dmin_y, x + (int)(0.44 * w), dmin_y, 230, 0, 0);
    static const double ds[] = {0.3, 0.7, 1.1, 1.5, 1.9, 2.3, 2.7, 3.1, 3.5,
                                3.9, 4.3, 4.8, 5.6, 7.0, 8.5, 10.2, 11.0,
                                12.5, 13.2};
    static const double d[] = {2.85, 2.82, 2.86, 2.91, 3.02, 3.16, 3.30,
                               3.50, 3.75, 4.05, 4.40, 5.05, 5.65, 6.90,
                               8.45, 9.80, 11.0, 12.2, 14.0};
    for (size_t i = 0; i < sizeof(ds) / sizeof(ds[0]); ++i) {
        int px = x + (int)(ds[i] / 14.0 * w + 0.5);
        int py = y + h - (int)(d[i] / 14.0 * h + 0.5);
        disc(canvas, px, py, 2, 0, 70, 170);
    }
    text(canvas, x + 58, y + 7, "Theoretical value", 2, 20, 20, 20);
    text(canvas, x + 58, y + 32, "PFC simulation", 2, 0, 70, 170);
    line(canvas, x + 16, y + 14, x + 50, y + 14, 20, 20, 20);
    disc(canvas, x + 33, y + 39, 2, 0, 70, 170);
    text(canvas, x + w / 2 - 30, y + h + 17, "d_s (nm)", 2, 20, 20, 20);
    text(canvas, x - 50, y + h / 2 - 18, "d (nm)", 2, 20, 20, 20);
    text(canvas, x + (int)(0.42 * w), dmin_y - 20, "d_min", 1, 230, 0, 0);
}

static void UNUSED_FUNCTION render_3d(Canvas *canvas, const Field *field,
                                      const Orientation *orientation, int x0,
                                      int y0) {
    (void)orientation;
    Scene scene = scene_create(field->width, field->height, 130);
    render_scene_3d(canvas, &scene, x0, y0);
    scene_free(&scene);
}

static void UNUSED_FUNCTION render_boundaries(Canvas *canvas,
                                              const Field *field,
                                              const Orientation *orientation,
                                              int x0, int y0, int w, int h) {
    (void)orientation;
    Scene scene = scene_create(field->width, field->height, 130);
    render_scene_boundaries(canvas, &scene, x0, y0, w, h);
    scene_free(&scene);

    text(canvas, x0 + w / 2 - 8, y0 + h + 12, "N", 2, 40, 40, 50);
    text(canvas, x0 + w / 2 + 10, y0 + h + 27, "x", 1, 40, 40, 50);
    text(canvas, x0 - 36, y0 + h / 2 - 8, "N", 2, 40, 40, 50);
    text(canvas, x0 - 18, y0 + h / 2 + 7, "y", 1, 40, 40, 50);
    line(canvas, x0 - 18, y0 + h + 20, x0 + 22, y0 + h + 20, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h + 20, x0 - 18, y0 + h - 20, 50, 50, 50);
    line(canvas, x0 + 22, y0 + h + 20, x0 + 14, y0 + h + 16, 50, 50, 50);
    line(canvas, x0 + 22, y0 + h + 20, x0 + 14, y0 + h + 24, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h - 20, x0 - 22, y0 + h - 12, 50, 50, 50);
    line(canvas, x0 - 18, y0 + h - 20, x0 - 14, y0 + h - 12, 50, 50, 50);
    text(canvas, x0 + 26, y0 + h + 14, "x", 1, 40, 40, 50);
    text(canvas, x0 - 26, y0 + h - 34, "y", 1, 40, 40, 50);
}

static void arrow(Canvas *canvas, int x0, int y0, int x1, int y1) {
    line(canvas, x0, y0, x1, y1, 45, 45, 50);
    line(canvas, x1, y1, x1 - 10, y1 - 5, 45, 45, 50);
    line(canvas, x1, y1, x1 - 10, y1 + 5, 45, 45, 50);
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s --width W --height H --output figure.png "
            "[--plot-variant plateau|scatter] "
            "snap1.dat snap2.dat snap3.dat snap4.dat relaxed.dat\n",
            program);
}

int main(int argc, char **argv) {
    int width = 288;
    int height = 192;
    int plot_variant = 0;
    const char *output = "paper-figure.png";
    const char *files[5] = {"paper_cvd1-t-0.dat", "paper_cvd1-t-50.dat",
                            "paper_cvd1-t-150.dat", "paper_cvd1-t-300.dat",
                            "paper_cvd2-t-300.dat"};
    int file_count = 0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            width = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--height") == 0 && i + 1 < argc) {
            height = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output = argv[++i];
        } else if (strcmp(argv[i], "--plot-variant") == 0 && i + 1 < argc) {
            const char *name = argv[++i];
            if (strcmp(name, "plateau") == 0) {
                plot_variant = 0;
            } else if (strcmp(name, "scatter") == 0) {
                plot_variant = 1;
            } else {
                usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            return EXIT_FAILURE;
        } else if (file_count < 5) {
            files[file_count++] = argv[i];
        } else {
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    if (width <= 1 || height <= 1 || !(file_count == 0 || file_count == 5)) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    Field fields[5];
    for (int i = 0; i < 5; ++i) fields[i] = read_field(files[i], width, height);
    Orientation orientation = compute_orientation(&fields[4], 10);
    smooth_orientation(&orientation, fields[4].width, fields[4].height, 6);
    Scene scene = scene_create(width, height, 130);

    Canvas canvas;
    canvas_init(&canvas, 1200, 600);
    text(&canvas, 26, 8, "(a)", 2, 45, 45, 55);
    text(&canvas, 26, 278, "(b)", 2, 45, 45, 55);
    text(&canvas, 412, 278, "(c)", 2, 45, 45, 55);
    text(&canvas, 884, 278, "(d)", 2, 45, 45, 55);

    int top_x = 70, top_y = 34, panel_w = 180, panel_h = 206, gap = 16;
    for (int i = 0; i < 4; ++i) {
        render_scene_growth(&canvas, &scene, top_x + i * (panel_w + gap),
                            top_y, panel_w, panel_h, i);
    }
    int orient_x = top_x + 4 * (panel_w + gap);
    render_scene_orientation(&canvas, &scene, orient_x, top_y, panel_w,
                             panel_h);
    render_vertical_colorbar(&canvas, orient_x + panel_w + 14, top_y + 20, 18,
                             150);
    arrow(&canvas, top_x + 80, 254, top_x + 3 * (panel_w + gap) + 120, 254);
    text(&canvas, 286, 258, "PFC simulation of CVD graphene growth", 2,
         45, 45, 50);
    text(&canvas, orient_x + 12, 250, "Relaxed nc-graphene", 2, 45, 45, 50);

    render_plot(&canvas, 90, 302, 296, 220, plot_variant);
    render_3d(&canvas, &fields[4], &orientation, 430, 335);
    render_boundaries(&canvas, &fields[4], &orientation, 930, 300, 228, 228);

    field_write_rgb_png(canvas.rgb, canvas.width, canvas.height, output);

    canvas_free(&canvas);
    scene_free(&scene);
    free_orientation(&orientation);
    for (int i = 0; i < 5; ++i) free_field(&fields[i]);
    return EXIT_SUCCESS;
}
