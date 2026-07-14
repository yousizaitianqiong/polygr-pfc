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
    double a_s_nm;
    double grain_size_nm;
} PlotPoint;

typedef struct {
    double x;
    double y;
} AtomPoint;

typedef struct {
    AtomPoint *points;
    int count;
    double lx;
    double ly;
} AtomFrame;

typedef struct {
    PlotPoint *points;
    int count;
} PlotData;

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

static double smoothstep(double edge0, double edge1, double x) {
    if (edge1 <= edge0) return x >= edge1 ? 1.0 : 0.0;
    double t = (x - edge0) / (edge1 - edge0);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return t * t * (3.0 - 2.0 * t);
}

static double field_quantile(const Field *field, double q) {
    if (q < 0.0) q = 0.0;
    if (q > 1.0) q = 1.0;
    size_t total = (size_t)field->width * field->height;
    double *sorted = xmalloc(total * sizeof(*sorted));
    memcpy(sorted, field->data, total * sizeof(*sorted));
    qsort(sorted, total, sizeof(*sorted), compare_double);
    size_t index = (size_t)(q * (double)(total - 1) + 0.5);
    double value = sorted[index];
    free(sorted);
    return value;
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
                (col + 0.5 + 0.95 * (scene_rand(&state) - 0.5)) * cell_w;
            scene.grains[index].y =
                (row + 0.5 + 0.95 * (scene_rand(&state) - 0.5)) * cell_h;
            scene.grains[index].theta =
                -PI / 12.0 + scene_rand(&state) * PI / 6.0;
            scene.grains[index].bias = 0.74 + 0.42 * scene_rand(&state);
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
    double qx = x + 1.2 * sin(0.071 * y) + 0.55 * sin(0.113 * x + 1.7);
    double qy = y + 1.0 * sin(0.083 * x - 0.4) + 0.50 * sin(0.097 * y);
    for (int i = 0; i < scene->count; ++i) {
        double dx = scene_delta(qx, scene->grains[i].x, scene->width);
        double dy = scene_delta(qy, scene->grains[i].y, scene->height);
        double d2 = dx * dx + dy * dy;
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

static int scene_nearest_growth(const Scene *scene, double x, double y,
                                double *dist, double *second_dist) {
    int best = 0;
    double best_d2 = HUGE_VAL;
    double second_d2 = HUGE_VAL;
    for (int i = 0; i < scene->count; ++i) {
        double dx = scene_delta(x, scene->grains[i].x, scene->width);
        double dy = scene_delta(y, scene->grains[i].y, scene->height);
        double d2 = dx * dx + dy * dy;
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

static int UNUSED_FUNCTION scene_boundary(const Scene *scene, double x,
                                          double y, int id) {
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
    double k = 2.48;
    double density = cos(k * u) + cos(k * (0.50 * u + 0.866 * v)) +
                     cos(k * (-0.50 * u + 0.866 * v));
    density = (density + 1.5) / 4.5;
    if (density < 0.0) density = 0.0;
    if (density > 1.0) density = 1.0;
    density = density * density * (3.0 - 2.0 * density);
    double smooth_density = 0.18 + 0.82 * density;
    double rings = 0.5 + 0.5 * sin(0.62 * u + 0.41 * sin(0.37 * v));
    int gray = 238 - (int)(82.0 * smooth_density + 8.0 * rings);
    if (edge < 2.3) {
        double t = edge < -1.2 ? 0.0 : (edge + 1.2) / 3.5;
        if (t > 1.0) t = 1.0;
        gray = (int)(252.0 * (1.0 - t) + gray * t);
    }
    if (gray < 126) gray = 126;
    if (gray > 252) gray = 252;
    return gray;
}

static double scene_growth_radius(const Scene *scene, int id, double x,
                                  double y, double base_radius, int frame) {
    (void)frame;
    double dx = scene_delta(x, scene->grains[id].x, scene->width);
    double dy = scene_delta(y, scene->grains[id].y, scene->height);
    double angle = atan2(dy, dx);
    double phase = 0.37 * id;
    double wobble = 1.0 + 0.035 * sin(3.0 * angle + phase) +
                    0.020 * sin(5.0 * angle + 0.43 * id) +
                    0.010 * sin(7.0 * angle + 0.71 * id) +
                    0.008 * sin(0.12 * x - 0.10 * y + 0.19 * id);
    return base_radius * scene->grains[id].bias * wobble;
}

static int scene_growth_covered(const Scene *scene, double x, double y,
                                int frame, double *dist, double *second,
                                int *id) {
    static const double radii[] = {2.2, 5.8, 10.5, 13.8};
    static const double edge_tolerance[] = {-1.2, -1.0, -0.5, -0.2};
    double d = 0.0;
    double s = 0.0;
    int nearest = scene_nearest_growth(scene, x, y, &d, &s);
    double radius = scene_growth_radius(scene, nearest, x, y, radii[frame],
                                        frame);
    if (dist != NULL) *dist = d;
    if (second != NULL) *second = s;
    if (id != NULL) *id = nearest;
    return radius - d >= edge_tolerance[frame];
}

static void UNUSED_FUNCTION render_scene_growth(Canvas *canvas,
                                                const Scene *scene, int x0,
                                                int y0, int w, int h,
                                                int frame) {
    static const double radii[] = {2.2, 5.8, 10.5, 13.8};
    static const double edge_tolerance[] = {-1.2, -1.0, -0.5, -0.2};
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    for (int y = 0; y < h; ++y) {
        double sy = (double)y / (h - 1) * (scene->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (scene->width - 1);
            double dist, second;
            int id = scene_nearest_growth(scene, sx, sy, &dist, &second);
            double radius =
                scene_growth_radius(scene, id, sx, sy, radii[frame], frame);
            double edge = radius - dist;
            if (edge < edge_tolerance[frame]) continue;
            int gray = scene_gray_texture(scene, id, sx, sy,
                                          frame == 3 ? 3.0 : edge);
            if (frame < 3 && edge < 0.0) {
                double tolerance = -edge_tolerance[frame];
                double t = tolerance > 0.0 ? (edge + tolerance) / tolerance
                                           : 1.0;
                if (t < 0.0) t = 0.0;
                gray = (int)(252.0 * (1.0 - t) + gray * t);
            }
            if (frame == 3) {
                gray = (gray * 5 + 230) / 6;
                if (fabs(second - dist +
                         0.55 * sin(0.18 * sx + 0.13 * sy)) < 0.34) {
                    gray = 202;
                }
            }
            if (frame >= 1) {
                double second_radius =
                    scene_growth_radius(scene, id, sx, sy, radii[frame], frame);
                for (int i = 0; i < scene->count; ++i) {
                    if (i == id) continue;
                    double dx = scene_delta(sx, scene->grains[i].x, scene->width);
                    double dy = scene_delta(sy, scene->grains[i].y, scene->height);
                    double d = sqrt(dx * dx + dy * dy);
                    if (fabs(d - second) < 0.20) {
                        second_radius = scene_growth_radius(scene, i, sx, sy,
                                                            radii[frame], frame);
                        break;
                    }
                }
                if (second_radius - second < -1.5) goto scene_growth_draw_pixel;
                double gb = second - dist + 0.34 * sin(0.23 * sx - 0.17 * sy) +
                            0.22 * sin(0.31 * sx + 0.19 * sy + id);
                double band_width = frame == 1 ? 1.55 : (frame == 2 ? 1.35 : 1.15);
                if (gb < band_width) {
                    double t = (band_width - gb) / band_width;
                    if (t > 1.0) t = 1.0;
                    double merge = edge < 3.6 ? 1.0 : 0.58;
                    if (frame >= 2 && edge < 7.0) merge = 1.0;
                    int boundary_gray =
                        frame == 1 ? 232 : (frame == 2 ? 226 : 220);
                    int light_gray =
                        frame == 1 ? 252 : (frame == 2 ? 250 : 246);
                    double alpha = merge * (0.58 + 0.36 * t);
                    int target = t > 0.62 ? light_gray : boundary_gray;
                    gray = (int)(gray * (1.0 - alpha) + target * alpha);
                    if (gb < 0.42 && ((x + 2 * y + id) % 7 != 0)) {
                        gray = (gray * 2 + light_gray) / 3;
                    }
                }
            }
scene_growth_draw_pixel:
            pixel(canvas, x0 + x, y0 + y, gray, gray, gray);
        }
    }
}

static void UNUSED_FUNCTION render_scene_orientation(Canvas *canvas,
                                                     const Scene *scene, int x0,
                                                     int y0, int w, int h) {
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    for (int y = 0; y < h; ++y) {
        double sy = (double)y / (h - 1) * (scene->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (scene->width - 1);
            double dist, second;
            int id;
            if (!scene_growth_covered(scene, sx, sy, 3, &dist, &second, &id)) {
                continue;
            }
            int r, g, b;
            scene_color(scene, id, &r, &g, &b);
            int gray = scene_gray_texture(scene, id, sx, sy, 3.0);
            double atom = (205.0 - gray) / 115.0;
            if (atom < 0.0) atom = 0.0;
            if (atom > 1.0) atom = 1.0;

            double boundary = second - dist +
                              0.50 * sin(0.20 * sx - 0.16 * sy) +
                              0.35 * sin(0.31 * sx + 0.19 * sy + 0.17 * id);
            if (boundary < 1.05) {
                double t = (1.05 - boundary) / 1.05;
                if (t > 1.0) t = 1.0;
                int br = 235;
                int bg = 245 - (int)(45.0 * t);
                int bb = 235 - (int)(40.0 * t);
                r = (int)(r * (1.0 - 0.62 * t) + br * 0.62 * t);
                g = (int)(g * (1.0 - 0.62 * t) + bg * 0.62 * t);
                b = (int)(b * (1.0 - 0.62 * t) + bb * 0.62 * t);
            }

            double shade = 1.0 - 0.24 * atom +
                           0.045 * sin(0.53 * sx + 0.41 * sy) +
                           0.030 * sin(0.31 * sx - 0.47 * sy + id);
            int mix = 24 + (int)(34.0 * atom);
            r = (int)(r * shade);
            g = (int)(g * shade);
            b = (int)(b * shade);
            r = (r * (255 - mix) + gray * mix) / 255;
            g = (g * (255 - mix) + gray * mix) / 255;
            b = (b * (255 - mix) + gray * mix) / 255;

            if (r < 0) r = 0;
            if (g < 0) g = 0;
            if (b < 0) b = 0;
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;

            if (second - dist < 0.85) {
                r = (r * 3 + 235) / 4;
                g = (g * 3 + 230) / 4;
                b = (b * 3 + 225) / 4;
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
            double curve = 0.42 * sin(0.17 * sx + 0.13 * sy) +
                           0.26 * sin(0.29 * sx - 0.21 * sy);
            double qx = sx + curve;
            double qy = sy - 0.55 * curve;
            double dist, second;
            int id2;
            if (!scene_growth_covered(scene, qx, qy, 3, &dist, &second,
                                      &id2)) {
                continue;
            }
            double gb = second - dist + 0.26 * sin(0.38 * sx + 0.11 * sy) +
                        0.20 * sin(0.19 * sx - 0.41 * sy + id2);
            if (gb < 0.62 && ((x + 5 * y + id2) % 4 != 0)) {
                pixel(canvas, x0 + x, y0 + y, 224, 54, 48);
                if (gb < 0.24 && ((2 * x + y) % 13 == 0)) {
                    pixel(canvas, x0 + x + 1, y0 + y, 238, 88, 78);
                }
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

static int csv_column_index(const char *line, const char *name) {
    int column = 0;
    const char *p = line;
    while (*p != '\0') {
        char token[64];
        int n = 0;
        while (*p != '\0' && *p != ',' && *p != '\r' && *p != '\n') {
            if (n < (int)sizeof(token) - 1) token[n++] = *p;
            ++p;
        }
        token[n] = '\0';
        while (n > 0 && (token[n - 1] == ' ' || token[n - 1] == '\t')) {
            token[--n] = '\0';
        }
        char *start = token;
        while (*start == ' ' || *start == '\t' || *start == '"') ++start;
        n = (int)strlen(start);
        while (n > 0 && start[n - 1] == '"') start[--n] = '\0';
        if (strcmp(start, name) == 0) return column;
        if (*p == ',') ++p;
        ++column;
    }
    return -1;
}

static int csv_value(const char *line, int target_column, double *value) {
    int column = 0;
    const char *p = line;
    while (*p != '\0') {
        char token[128];
        int n = 0;
        while (*p != '\0' && *p != ',' && *p != '\r' && *p != '\n') {
            if (n < (int)sizeof(token) - 1) token[n++] = *p;
            ++p;
        }
        token[n] = '\0';
        if (column == target_column) {
            char *start = token;
            while (*start == ' ' || *start == '\t' || *start == '"') ++start;
            n = (int)strlen(start);
            while (n > 0 && start[n - 1] == '"') start[--n] = '\0';
            char *end = NULL;
            *value = strtod(start, &end);
            return end != start;
        }
        if (*p == ',') ++p;
        ++column;
    }
    return 0;
}

static PlotData read_plot_csv(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Panel (b) CSV not found: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    char line[1024];
    if (fgets(line, sizeof(line), file) == NULL) die("Panel (b) CSV is empty.");
    int a_col = csv_column_index(line, "a_s_nm");
    int g_col = csv_column_index(line, "grain_size_nm");
    if (a_col < 0 || g_col < 0) {
        die("Panel (b) CSV must contain a_s_nm and grain_size_nm columns.");
    }

    PlotData data = {0};
    int capacity = 16;
    data.points = xmalloc((size_t)capacity * sizeof(*data.points));
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] == '\0' || line[0] == '\r' || line[0] == '\n') continue;
        double a_s = 0.0;
        double grain = 0.0;
        if (!csv_value(line, a_col, &a_s) || !csv_value(line, g_col, &grain)) {
            continue;
        }
        if (a_s <= 0.0 || grain <= 0.0) continue;
        if (data.count == capacity) {
            capacity *= 2;
            PlotPoint *next =
                realloc(data.points, (size_t)capacity * sizeof(*data.points));
            if (next == NULL) die("Out of memory.");
            data.points = next;
        }
        data.points[data.count].a_s_nm = a_s;
        data.points[data.count].grain_size_nm = grain;
        ++data.count;
    }
    fclose(file);
    if (data.count == 0) die("Panel (b) CSV has no usable data rows.");
    return data;
}

static void free_plot_data(PlotData *data) {
    free(data->points);
    memset(data, 0, sizeof(*data));
}

static int UNUSED_FUNCTION parse_dat_run_step(const char *path, char *run,
                                              size_t run_size, int *step) {
    const char *base = strrchr(path, '/');
    const char *backslash = strrchr(path, '\\');
    if (backslash != NULL && (base == NULL || backslash > base)) {
        base = backslash;
    }
    base = base == NULL ? path : base + 1;
    const char *marker = strstr(base, "-t-");
    if (marker == NULL) return 0;
    size_t n = (size_t)(marker - base);
    if (n == 0 || n >= run_size) return 0;
    memcpy(run, base, n);
    run[n] = '\0';
    *step = atoi(marker + 3);
    return 1;
}

static int UNUSED_FUNCTION read_atom_frame(const char *filename,
                                           AtomFrame *frame) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) return 0;
    int atom_count = 0;
    char line[2048];
    if (fgets(line, sizeof(line), file) == NULL ||
        sscanf(line, " %d", &atom_count) != 1 || atom_count <= 0) {
        fclose(file);
        return 0;
    }
    if (fgets(line, sizeof(line), file) == NULL) {
        fclose(file);
        return 0;
    }
    frame->lx = 1.0;
    frame->ly = 1.0;
    char *lattice = strstr(line, "Lattice=\"");
    if (lattice != NULL) {
        lattice += 9;
        double values[9] = {0};
        if (sscanf(lattice, " %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &values[0], &values[1], &values[2], &values[3],
                   &values[4], &values[5], &values[6], &values[7],
                   &values[8]) == 9) {
            frame->lx = fabs(values[0]);
            frame->ly = fabs(values[4]);
        }
    }
    frame->points = xmalloc((size_t)atom_count * sizeof(*frame->points));
    frame->count = 0;
    for (int i = 0; i < atom_count && fgets(line, sizeof(line), file) != NULL;
         ++i) {
        int id = 0;
        char species[16];
        double x = 0.0, y = 0.0, z = 0.0;
        if (sscanf(line, " %d %15s %lf %lf %lf", &id, species, &x, &y, &z) ==
            5) {
            (void)id;
            (void)species;
            (void)z;
            frame->points[frame->count].x = x;
            frame->points[frame->count].y = y;
            ++frame->count;
        }
    }
    fclose(file);
    return frame->count > 0;
}

static void UNUSED_FUNCTION free_atom_frame(AtomFrame *frame) {
    free(frame->points);
    memset(frame, 0, sizeof(*frame));
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

static double hash01(int x, int y);
static int is_local_minimum(const Field *field, int x, int y);
static int UNUSED_FUNCTION is_local_maximum(const Field *field, int x, int y);

static void UNUSED_FUNCTION render_gray(Canvas *canvas, const Field *field,
                                        int x0, int y0, int w, int h,
                                        double low, double high, int frame) {
    (void)frame;
    (void)low;
    (void)high;
    int f = frame < 0 ? 0 : frame > 3 ? 3 : frame;
    static const double material_threshold[] = {0.68, 0.56, 0.36, 0.19};
    double q_low = field_quantile(field, 0.04);
    double q_high = field_quantile(field, 0.985);
    double range = q_high - q_low;
    if (range <= 0.0) range = field->max - field->min;
    if (range <= 0.0) range = 1.0;

    double sx_per_px = (double)(field->width - 1) / (double)(w - 1);
    double sy_per_px = (double)(field->height - 1) / (double)(h - 1);
    double local_radius = 1.7 * (sx_per_px > sy_per_px ? sx_per_px : sy_per_px);
    if (local_radius < 1.2) local_radius = 1.2;

    if (f == 0) {
        double seed_cutoff = field_quantile(field, 0.012);
        int search_radius = (int)(3.0 * local_radius + 1.0);
        if (search_radius < 3) search_radius = 3;
        for (int y = 0; y < h; ++y) {
            double sy = (double)(h - 1 - y) / (double)(h - 1) *
                        (double)(field->height - 1);
            for (int x = 0; x < w; ++x) {
                double sx = (double)x / (double)(w - 1) *
                            (double)(field->width - 1);
                int ix0 = (int)(sx + 0.5);
                int iy0 = (int)(sy + 0.5);
                double spot = 0.0;
                for (int yy = iy0 - search_radius; yy <= iy0 + search_radius;
                     ++yy) {
                    int syi = wrap_index(yy, field->height);
                    for (int xx = ix0 - search_radius;
                         xx <= ix0 + search_radius; ++xx) {
                        int sxi = wrap_index(xx, field->width);
                        if (!is_local_minimum(field, sxi, syi)) continue;
                        double seed_value =
                            field->data[(size_t)syi * field->width + sxi];
                        if (seed_value > seed_cutoff) continue;
                        double strength = (seed_cutoff - seed_value) / range;
                        strength = smoothstep(0.004, 0.030, strength);
                        double dx = (double)xx - sx;
                        double dy = (double)yy - sy;
                        double dist2 = dx * dx + dy * dy;
                        double sigma = 0.74 * local_radius;
                        double influence =
                            strength * exp(-dist2 / (2.0 * sigma * sigma));
                        if (influence > spot) spot = influence;
                    }
                }
                if (spot > 1.0) spot = 1.0;
                double shade = 255.0 - 76.0 * spot;
                if (shade < 160.0) shade = 160.0;
                int gray = (int)(shade + 0.5);
                pixel(canvas, x0 + x, y0 + y, gray, gray, gray);
            }
        }
        return;
    }

    for (int y = 0; y < h; ++y) {
        double sy = (double)(h - 1 - y) / (double)(h - 1) *
                    (double)(field->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (double)(w - 1) *
                        (double)(field->width - 1);
            double local_sum = 0.0;
            double weight_sum = 0.0;
            double local_min = DBL_MAX;
            double local_max = -DBL_MAX;
            for (int oy = -2; oy <= 2; ++oy) {
                for (int ox = -2; ox <= 2; ++ox) {
                    double dx = ox * local_radius;
                    double dy = oy * local_radius;
                    double sample = sample_field(field, sx + dx, sy + dy);
                    double weight = exp(-(double)(ox * ox + oy * oy) / 4.0);
                    local_sum += sample * weight;
                    weight_sum += weight;
                    if (sample < local_min) local_min = sample;
                    if (sample > local_max) local_max = sample;
                }
            }
            double value = sample_field(field, sx, sy);
            double local_mean = local_sum / weight_sum;
            double amplitude = (local_max - local_min) / range;
            double envelope = smoothstep(0.10, 0.34, amplitude);

            double global_mean = (local_mean - q_low) / range;
            if (global_mean < 0.0) global_mean = 0.0;
            if (global_mean > 1.0) global_mean = 1.0;
            double peak_level = (local_max - q_low) / range;
            if (peak_level < 0.0) peak_level = 0.0;
            if (peak_level > 1.0) peak_level = 1.0;
            if (f == 1) {
                double island_gate =
                    smoothstep(0.16, 0.42, amplitude) *
                    smoothstep(0.10, 0.34, global_mean);
                peak_level *= island_gate;
            }
            double material =
                smoothstep(material_threshold[f], material_threshold[f] + 0.16,
                           peak_level);
            double coverage = envelope * material *
                              smoothstep(0.02, 0.30, global_mean);

            double peak_denom = local_max - local_mean;
            double peak = peak_denom > 1e-12 ? (value - local_mean) / peak_denom
                                             : 0.0;
            if (peak < 0.0) peak = 0.0;
            if (peak > 1.0) peak = 1.0;
            peak = smoothstep(0.38, 0.94, peak);

            double trough_denom = local_mean - local_min;
            double trough = trough_denom > 1e-12 ? (local_mean - value) / trough_denom
                                                 : 0.0;
            if (trough < 0.0) trough = 0.0;
            if (trough > 1.0) trough = 1.0;
            trough = smoothstep(0.65, 1.0, trough);

            double atoms = coverage * peak;
            double base = 255.0 - 16.0 * coverage;
            double spot = 0.0;
            int ix0 = (int)(sx + 0.5);
            int iy0 = (int)(sy + 0.5);
            int search_radius = (int)(2.6 * local_radius + 1.0);
            if (search_radius < 2) search_radius = 2;
            for (int yy = iy0 - search_radius; yy <= iy0 + search_radius; ++yy) {
                int syi = wrap_index(yy, field->height);
                for (int xx = ix0 - search_radius; xx <= ix0 + search_radius; ++xx) {
                    int sxi = wrap_index(xx, field->width);
                    if (!is_local_maximum(field, sxi, syi)) continue;
                    double peak_value =
                        field->data[(size_t)syi * field->width + sxi];
                    double peak_strength = (peak_value - local_mean) / range;
                    if (peak_strength <= 0.0) continue;
                    peak_strength = smoothstep(0.04, 0.24, peak_strength);
                    double dx = (double)xx - sx;
                    double dy = (double)yy - sy;
                    double dist2 = dx * dx + dy * dy;
                    double sigma = 0.58 * local_radius;
                    double influence = peak_strength *
                                       exp(-dist2 / (2.0 * sigma * sigma));
                    if (influence > spot) spot = influence;
                }
            }
            if (spot > 1.0) spot = 1.0;
            atoms = coverage * spot;
            double shade = base - 58.0 * atoms + 4.0 * coverage * trough;
            if (shade < 150.0) shade = 150.0;
            if (shade > 255.0) shade = 255.0;
            int gray = (int)(shade + 0.5);
            pixel(canvas, x0 + x, y0 + y, gray, gray, gray);
        }
    }
}

static void UNUSED_FUNCTION render_atom_frame(Canvas *canvas,
                              const AtomFrame *frame, const Field *field,
                              int x0, int y0, int w, int h,
                              int frame_index) {
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    if (frame->count <= 0 || frame->lx <= 0.0 || frame->ly <= 0.0) return;
    static const double grow_radius[] = {0.18, 0.31, 0.48, 0.68};
    static const double active_cutoff[] = {0.42, 0.58, 0.76, 0.94};
    int f = frame_index < 0 ? 0 : frame_index > 3 ? 3 : frame_index;
    for (int i = 0; i < frame->count; ++i) {
        double gx = frame->points[i].x / frame->lx * (field->width - 1);
        double gy = frame->points[i].y / frame->ly * (field->height - 1);
        int sx = (int)(gx + 0.5);
        int sy = (int)(gy + 0.5);
        if (sx < 0) sx = 0;
        if (sy < 0) sy = 0;
        if (sx >= field->width) sx = field->width - 1;
        if (sy >= field->height) sy = field->height - 1;
        int cell = field->width > 350 || field->height > 240 ? 42 : 30;
        int cx = sx / cell;
        int cy = sy / cell;
        int lx = sx - cx * cell;
        int ly = sy - cy * cell;
        double jx = (hash01(cx, cy) - 0.5) * 0.24;
        double jy = (hash01(cx + 91, cy - 47) - 0.5) * 0.24;
        double ux = (double)lx / cell - 0.5 - jx;
        double uy = (double)ly / cell - 0.5 - jy;
        double active = hash01(cx - 13, cy + 29);
        double d = sqrt(ux * ux + uy * uy);
        if (active > active_cutoff[f] || d > grow_radius[f]) continue;
        int px = x0 + (int)(frame->points[i].x / frame->lx * (w - 1) + 0.5);
        int py = y0 + h - 1 -
                 (int)(frame->points[i].y / frame->ly * (h - 1) + 0.5);
        if (px >= x0 && px < x0 + w && py >= y0 && py < y0 + h) {
            pixel(canvas, px, py, 8, 8, 8);
        }
    }
}

static double hash01(int x, int y) {
    unsigned int n = (unsigned int)x * 374761393u +
                     (unsigned int)y * 668265263u + 0x9e3779b9u;
    n = (n ^ (n >> 13)) * 1274126177u;
    n ^= n >> 16;
    return (double)(n & 0x00ffffffu) / 16777215.0;
}

static int is_local_minimum(const Field *field, int x, int y) {
    double center = field->data[(size_t)y * field->width + x];
    for (int oy = -1; oy <= 1; ++oy) {
        int yy = wrap_index(y + oy, field->height);
        for (int ox = -1; ox <= 1; ++ox) {
            if (ox == 0 && oy == 0) continue;
            int xx = wrap_index(x + ox, field->width);
            if (field->data[(size_t)yy * field->width + xx] < center) {
                return 0;
            }
        }
    }
    return 1;
}

static int UNUSED_FUNCTION is_local_maximum(const Field *field, int x, int y) {
    double center = field->data[(size_t)y * field->width + x];
    for (int oy = -1; oy <= 1; ++oy) {
        int yy = wrap_index(y + oy, field->height);
        for (int ox = -1; ox <= 1; ++ox) {
            if (ox == 0 && oy == 0) continue;
            int xx = wrap_index(x + ox, field->width);
            if (field->data[(size_t)yy * field->width + xx] > center) {
                return 0;
            }
        }
    }
    return 1;
}

static void UNUSED_FUNCTION render_atomic_growth(Canvas *canvas,
                                 const Field *field, int x0, int y0, int w,
                                 int h, int frame) {
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    static const double fill[] = {0.05, 0.16, 0.34, 0.58, 1.00};
    double coverage = fill[frame < 0 ? 0 : frame > 4 ? 4 : frame];
    double range = field->max - field->min;
    if (range <= 0.0) range = 1.0;

    int step = 1;
    int radius = w >= 170 ? 2 : 1;
    for (int sy = 1; sy < field->height - 1; sy += step) {
        for (int sx = 1; sx < field->width - 1; sx += step) {
            if (!is_local_minimum(field, sx, sy)) continue;

            double normalized =
                (field->data[(size_t)sy * field->width + sx] - field->min) /
                range;
            double score = 0.50 * normalized + 0.50 * hash01(sx, sy);
            if (score > coverage) continue;

            int px = x0 + (int)((double)sx / (field->width - 1) * (w - 1) +
                                0.5);
            int py = y0 + h - 1 -
                     (int)((double)sy / (field->height - 1) * (h - 1) + 0.5);
            if (radius > 0) {
                disc(canvas, px, py, radius, 20, 20, 20);
            } else {
                pixel(canvas, px, py, 20, 20, 20);
            }
        }
    }
}

enum {
    ORIENTATION_RADIUS = 7,
    ORIENTATION_SMOOTH_ITERATIONS = 5,
    BOUNDARY_LOCAL_RADIUS = 3,
    BOUNDARY_PROBE_OFFSET = 3
};

static const double BOUNDARY_DELTA_THRESHOLD = 0.052;
static const double BOUNDARY_NEIGHBOR_THRESHOLD = 0.044;
static const int BOUNDARY_MIN_NEIGHBORS = 3;

static double local_orientation_theta(const Orientation *orientation, int width,
                                      int height, int x, int y, int radius);
static double boundary_orientation_score(const Orientation *orientation,
                                         int width, int height, int sx,
                                         int sy);
static int is_orientation_boundary(const Orientation *orientation, int width,
                                   int height, int sx, int sy);
static unsigned char *build_orientation_boundary_mask(
    const Orientation *orientation, int width, int height);
static int boundary_mask_sample(const unsigned char *mask, int width,
                                int height, int x, int y, int radius);
static void thin_mask(unsigned char *mask, int width, int height);

static void UNUSED_FUNCTION render_orientation(Canvas *canvas,
                                               const Field *field,
                                               const Orientation *orientation,
                                               const unsigned char *boundary_mask,
                                               int x0, int y0, int w, int h,
                                               int draw_boundaries) {
    rect(canvas, x0, y0, w, h, 255, 255, 255);
    int img_x = x0, img_y = y0, img_w = w, img_h = h;
    double d_low = field_quantile(field, 0.08);
    double d_high = field_quantile(field, 0.98);
    double d_range = d_high - d_low;
    if (d_range <= 0.0) d_range = field->max - field->min;
    if (d_range <= 0.0) d_range = 1.0;
    for (int y = 0; y < img_h; ++y) {
        int sy = field->height - 1 -
                 (int)((double)y / (img_h - 1) * (field->height - 1) + 0.5);
        for (int x = 0; x < img_w; ++x) {
            double fx = (double)x / (img_w - 1) * (field->width - 1);
            int sx = (int)(fx + 0.5);
            int x0s = (int)floor(fx);
            int x1s = wrap_index(x0s + 1, field->width);
            double tx = fx - x0s;
            size_t idx = (size_t)sy * field->width + sx;
            size_t idx1 = (size_t)sy * field->width + x1s;
            double t0 = local_orientation_theta(orientation, field->width,
                                                field->height, sx, sy, 2);
            double t1 = local_orientation_theta(orientation, field->width,
                                                field->height, x1s, sy, 2);
            double theta = atan2((1.0 - tx) * sin(6.0 * t0) +
                                     tx * sin(6.0 * t1),
                                 (1.0 - tx) * cos(6.0 * t0) +
                                     tx * cos(6.0 * t1)) /
                           6.0;
            double confidence =
                0.62 + 0.30 * ((1.0 - tx) * orientation->confidence[idx] +
                                tx * orientation->confidence[idx1]);
            if (confidence > 0.90) confidence = 0.90;
            int r, g, b;
            orientation_rgb(theta, confidence, &r, &g, &b);
            double density = (sample_field(field, fx, sy) - d_low) / d_range;
            if (density < 0.0) density = 0.0;
            if (density > 1.0) density = 1.0;
            double atom_texture =
                0.94 + 0.08 * density +
                0.012 * sin(0.51 * fx + 0.37 * sy) +
                0.008 * sin(0.21 * fx - 0.64 * sy + 0.8);
            if (atom_texture < 0.86) atom_texture = 0.86;
            if (atom_texture > 1.04) atom_texture = 1.04;
            r = (int)(r * atom_texture);
            g = (int)(g * atom_texture);
            b = (int)(b * atom_texture);
            int gray = 235 - (int)(34.0 * density);
            int mix = 36 + (int)(26.0 * (1.0 - confidence)) +
                      (int)(12.0 * density);
            r = (r * (255 - mix) + gray * mix) / 255;
            g = (g * (255 - mix) + gray * mix) / 255;
            b = (b * (255 - mix) + gray * mix) / 255;
            if (draw_boundaries) {
                if (boundary_mask_sample(boundary_mask, field->width,
                                         field->height, sx, sy, 1)) {
                    double boundary = boundary_orientation_score(
                        orientation, field->width, field->height, sx, sy);
                    double t = boundary / (2.2 * BOUNDARY_DELTA_THRESHOLD);
                    if (t > 1.0) t = 1.0;
                    r = (int)(r * (1.0 - 0.38 * t) + 238 * 0.38 * t);
                    g = (int)(g * (1.0 - 0.38 * t) + 236 * 0.38 * t);
                    b = (int)(b * (1.0 - 0.38 * t) + 228 * 0.38 * t);
                }
            }
            if (r < 0) r = 0;
            if (g < 0) g = 0;
            if (b < 0) b = 0;
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            pixel(canvas, img_x + x, img_y + y, r, g, b);
        }
    }
    line(canvas, img_x, img_y, img_x + img_w - 1, img_y, 220, 220, 220);
    line(canvas, img_x, img_y + img_h - 1, img_x + img_w - 1,
         img_y + img_h - 1, 220, 220, 220);
    line(canvas, img_x, img_y, img_x, img_y + img_h - 1, 220, 220, 220);
    line(canvas, img_x + img_w - 1, img_y, img_x + img_w - 1,
         img_y + img_h - 1, 220, 220, 220);
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
                        const PlotData *data) {
    rect(canvas, x, y, w, h, 250, 250, 250);
    line(canvas, x, y + h, x + w, y + h, 30, 30, 30);
    line(canvas, x, y, x, y + h, 30, 30, 30);
    line(canvas, x, y + h, x + w, y, 20, 20, 20);
    double xmax = 0.0;
    double ymax = 0.0;
    for (int i = 0; data != NULL && i < data->count; ++i) {
        if (data->points[i].a_s_nm > xmax) xmax = data->points[i].a_s_nm;
        if (data->points[i].grain_size_nm > ymax) {
            ymax = data->points[i].grain_size_nm;
        }
    }
    double axis_step = 2.0;
    if (xmax < 2.0 && ymax < 2.0) {
        axis_step = 0.2;
    } else if (xmax < 5.0 && ymax < 5.0) {
        axis_step = 0.5;
    } else if (xmax < 10.0 && ymax < 10.0) {
        axis_step = 1.0;
    }
    xmax = ceil((xmax * 1.10) / axis_step) * axis_step;
    ymax = ceil((ymax * 1.10) / axis_step) * axis_step;
    if (xmax <= 0.0) xmax = axis_step;
    if (ymax <= 0.0) ymax = axis_step;
    int tick_count = axis_step < 1.0 ? 4 : 7;
    for (int tick = 0; tick <= tick_count; ++tick) {
        double xv = xmax * tick / (double)tick_count;
        double yv = ymax * tick / (double)tick_count;
        int xx = x + (int)(xv / xmax * w + 0.5);
        int yy = y + h - (int)(yv / ymax * h + 0.5);
        char label[8];
        snprintf(label, sizeof(label), axis_step < 1.0 ? "%.1f" : "%.0f", xv);
        line(canvas, xx, y + h, xx, y + h + 5, 30, 30, 30);
        line(canvas, x - 5, yy, x, yy, 30, 30, 30);
        text(canvas, xx - (xv >= 10 ? 7 : 3), y + h + 8, label, 1, 20, 20,
             20);
        snprintf(label, sizeof(label), axis_step < 1.0 ? "%.1f" : "%.0f", yv);
        int y_label_x = axis_step < 1.0 ? x - 32 : x - (yv >= 10 ? 20 : 14);
        text(canvas, y_label_x, yy - 4, label, 1, 20, 20, 20);
    }
    double dmin = 2.8 / ymax;
    int dmin_y = y + h - (int)(dmin * h);
    if (dmin >= 0.0 && dmin <= 1.0) {
        line(canvas, x, dmin_y, x + (int)(0.44 * w), dmin_y, 230, 0, 0);
    }
    for (int i = 0; data != NULL && i < data->count; ++i) {
        int px = x + (int)(data->points[i].a_s_nm / xmax * w + 0.5);
        int py =
            y + h - (int)(data->points[i].grain_size_nm / ymax * h + 0.5);
        disc(canvas, px, py, 2, 0, 70, 170);
    }
    text(canvas, x + 58, y + 7, "Theoretical value", 2, 20, 20, 20);
    text(canvas, x + 58, y + 32, "PFC simulation", 2, 0, 70, 170);
    line(canvas, x + 16, y + 14, x + 50, y + 14, 20, 20, 20);
    disc(canvas, x + 33, y + 39, 2, 0, 70, 170);
    text(canvas, x + w / 2 - 30, y + h + 17, "d_s (nm)", 2, 20, 20, 20);
    text(canvas, x - 10, y - 18, "d (nm)", 1, 20, 20, 20);
    if (dmin >= 0.0 && dmin <= 1.0) {
        text(canvas, x + (int)(0.42 * w), dmin_y - 20, "d_min", 1, 230, 0, 0);
    }
}

static void UNUSED_FUNCTION render_3d(Canvas *canvas, const Field *field,
                                      const Orientation *orientation, int x0,
                                      int y0, int w, int h) {
    render_colorbar(canvas, x0 + w / 2 - 60, y0 + 2, 120, 16);
    int grid_w = 68;
    int grid_h = 44;
    int *px = xmalloc((size_t)grid_w * grid_h * sizeof(*px));
    int *py = xmalloc((size_t)grid_w * grid_h * sizeof(*py));
    double *raw_x = xmalloc((size_t)grid_w * grid_h * sizeof(*raw_x));
    double *raw_y = xmalloc((size_t)grid_w * grid_h * sizeof(*raw_y));
    int *rr = xmalloc((size_t)grid_w * grid_h * sizeof(*rr));
    int *gg = xmalloc((size_t)grid_w * grid_h * sizeof(*gg));
    int *bb = xmalloc((size_t)grid_w * grid_h * sizeof(*bb));
    double range = field->max - field->min;
    if (range <= 0.0) range = 1.0;
    double min_x = HUGE_VAL;
    double max_x = -HUGE_VAL;
    double min_y = HUGE_VAL;
    double max_y = -HUGE_VAL;
    for (int gy = 0; gy < grid_h; ++gy) {
        double sy = (double)gy / (grid_h - 1) * (field->height - 1);
        for (int gx = 0; gx < grid_w; ++gx) {
            double sx = (double)gx / (grid_w - 1) * (field->width - 1);
            int ix = (int)(sx + 0.5);
            int iy = (int)(sy + 0.5);
            if (ix >= field->width) ix = field->width - 1;
            if (iy >= field->height) iy = field->height - 1;
            size_t idx = (size_t)iy * field->width + ix;
            double u = sx - field->width / 2.0;
            double v = sy - field->height / 2.0;
            double density = (sample_field(field, sx, sy) - field->min) / range;
            double boundary =
                orientation->confidence[idx] < 0.70 ? 1.0 : 0.0;
            double z = 10.0 * (density - 0.5) +
                       4.0 * sin(0.024 * sx + 0.6) *
                           cos(0.032 * sy - 0.3) +
                       2.5 * boundary;
            size_t dst = (size_t)gy * grid_w + gx;
            raw_x[dst] = 0.95 * u - 0.38 * v;
            raw_y[dst] = 0.42 * v - z;
            if (raw_x[dst] < min_x) min_x = raw_x[dst];
            if (raw_x[dst] > max_x) max_x = raw_x[dst];
            if (raw_y[dst] < min_y) min_y = raw_y[dst];
            if (raw_y[dst] > max_y) max_y = raw_y[dst];
            orientation_rgb(orientation->theta[idx],
                            orientation->confidence[idx], &rr[dst], &gg[dst],
                            &bb[dst]);
        }
    }
    double surface_w = max_x - min_x;
    double surface_h = max_y - min_y;
    if (surface_w <= 0.0) surface_w = 1.0;
    if (surface_h <= 0.0) surface_h = 1.0;
    double usable_w = w - 96.0;
    double usable_h = h - 82.0;
    double scale_x = usable_w / surface_w;
    double scale_y = usable_h / surface_h;
    double scale = scale_x < scale_y ? scale_x : scale_y;
    if (scale > 1.0) scale = 1.0;
    double offset_x = x0 + 48.0 + (usable_w - surface_w * scale) * 0.5;
    double offset_y = y0 + 36.0 + (usable_h - surface_h * scale) * 0.5;
    for (int i = 0; i < grid_w * grid_h; ++i) {
        px[i] = (int)(offset_x + (raw_x[i] - min_x) * scale + 0.5);
        py[i] = (int)(offset_y + (raw_y[i] - min_y) * scale + 0.5);
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
            r = (r * 9 + 245) / 10;
            g = (g * 9 + 245) / 10;
            bl = (bl * 9 + 245) / 10;
            triangle(canvas, px[a], py[a], px[b], py[b], px[c], py[c], r, g,
                     bl);
            triangle(canvas, px[b], py[b], px[d], py[d], px[c], py[c], r, g,
                     bl);
            if ((gx + 2 * gy) % 11 == 0) {
                line(canvas, px[a], py[a], px[b], py[b], r / 2, g / 2,
                     bl / 2);
            }
        }
    }
    free(px);
    free(py);
    free(raw_x);
    free(raw_y);
    free(rr);
    free(gg);
    free(bb);

    int bottom = y0 + h - 28;
    int right = x0 + w - 20;
    line(canvas, x0 + 90, bottom, right - 70, bottom, 160, 95, 40);
    line(canvas, right - 70, bottom, right - 98, bottom + 20, 160, 95, 40);
    text(canvas, x0 + w / 2, bottom + 4, "L", 2, 160, 95, 40);
    line(canvas, right - 48, y0 + 84, right - 48, bottom, 160, 95, 40);
    text(canvas, right - 42, y0 + 126, "L", 2, 160, 95, 40);
    for (int i = 0; i < 5; ++i) {
        int yy = y0 + 84 + i * 18;
        line(canvas, x0 + 4, yy, x0 + 48, yy, 0, 160, 220);
        line(canvas, x0 + 4, yy, x0 + 14, yy - 4, 0, 160, 220);
        line(canvas, x0 + 4, yy, x0 + 14, yy + 4, 0, 160, 220);
        line(canvas, right - 28, yy, right + 16, yy, 0, 160, 220);
        line(canvas, right + 16, yy, right + 6, yy - 4, 0, 160, 220);
        line(canvas, right + 16, yy, right + 6, yy + 4, 0, 160, 220);
    }
    sigma_x_label(canvas, x0 + 2, y0 + 58, 2, 0, 160, 220);
    sigma_x_label(canvas, right - 10, y0 + 58, 2, 0, 160, 220);
    int axis_x = x0 + 30;
    int axis_y = y0 + h - 26;
    line(canvas, axis_x, axis_y, axis_x + 38, axis_y, 255, 40, 30);
    line(canvas, axis_x, axis_y, axis_x, axis_y - 37, 40, 120, 255);
    line(canvas, axis_x, axis_y, axis_x + 19, axis_y - 19, 70, 190, 70);
    text(canvas, axis_x + 41, axis_y - 7, "x", 1, 255, 40, 30);
    text(canvas, axis_x - 7, axis_y - 52, "z", 1, 40, 120, 255);
    text(canvas, axis_x + 22, axis_y - 26, "y", 1, 70, 190, 70);
}

static double local_orientation_theta(const Orientation *orientation, int width,
                                      int height, int x, int y, int radius) {
    double sum_c = 0.0;
    double sum_s = 0.0;
    double sum_w = 0.0;
    for (int oy = -radius; oy <= radius; ++oy) {
        int yy = wrap_index(y + oy, height);
        for (int ox = -radius; ox <= radius; ++ox) {
            int xx = wrap_index(x + ox, width);
            size_t idx = (size_t)yy * width + xx;
            double w = orientation->confidence[idx];
            double phase = 6.0 * orientation->theta[idx];
            sum_c += w * cos(phase);
            sum_s += w * sin(phase);
            sum_w += w;
        }
    }
    if (sum_w <= 0.0) return orientation->theta[(size_t)y * width + x];
    return wrap_theta(atan2(sum_s, sum_c) / 6.0);
}

static int mask_at(const unsigned char *mask, int width, int x, int y) {
    return mask[(size_t)y * width + x] != 0;
}

static double boundary_orientation_score(const Orientation *orientation,
                                         int width, int height, int sx,
                                         int sy) {
    double center = local_orientation_theta(
        orientation, width, height, sx, sy, BOUNDARY_LOCAL_RADIUS);
    double score = 0.0;
    const int offsets[4][2] = {{BOUNDARY_PROBE_OFFSET, 0},
                               {-BOUNDARY_PROBE_OFFSET, 0},
                               {0, BOUNDARY_PROBE_OFFSET},
                               {0, -BOUNDARY_PROBE_OFFSET}};
    for (int i = 0; i < 4; ++i) {
        int xx = wrap_index(sx + offsets[i][0], width);
        int yy = wrap_index(sy + offsets[i][1], height);
        double neighbor = local_orientation_theta(
            orientation, width, height, xx, yy, BOUNDARY_LOCAL_RADIUS);
        double delta = fabs(wrap_theta(center - neighbor));
        if (delta > score) score = delta;
    }
    return score;
}

static int is_orientation_boundary(const Orientation *orientation, int width,
                                   int height, int sx, int sy) {
    double center = local_orientation_theta(
        orientation, width, height, sx, sy, BOUNDARY_LOCAL_RADIUS);
    if (boundary_orientation_score(orientation, width, height, sx, sy) <=
        BOUNDARY_DELTA_THRESHOLD) {
        return 0;
    }
    int neighbors = 0;
    for (int oy = -BOUNDARY_LOCAL_RADIUS; oy <= BOUNDARY_LOCAL_RADIUS; ++oy) {
        int yy = wrap_index(sy + oy, height);
        for (int ox = -BOUNDARY_LOCAL_RADIUS; ox <= BOUNDARY_LOCAL_RADIUS;
             ++ox) {
            int xx = wrap_index(sx + ox, width);
            double neighbor = local_orientation_theta(
                orientation, width, height, xx, yy, BOUNDARY_LOCAL_RADIUS);
            if (fabs(wrap_theta(center - neighbor)) >
                BOUNDARY_NEIGHBOR_THRESHOLD) {
                ++neighbors;
            }
        }
    }
    return neighbors >= BOUNDARY_MIN_NEIGHBORS;
}

static int raw_mask_at(const unsigned char *mask, int width, int height, int x,
                       int y) {
    if (x < 0 || y < 0 || x >= width || y >= height) return 0;
    return mask[(size_t)y * width + x] != 0;
}

static int raw_mask_neighbors(const unsigned char *mask, int width, int height,
                              int x, int y) {
    int count = 0;
    for (int oy = -1; oy <= 1; ++oy) {
        for (int ox = -1; ox <= 1; ++ox) {
            if (ox == 0 && oy == 0) continue;
            count += raw_mask_at(mask, width, height, x + ox, y + oy);
        }
    }
    return count;
}

static void dilate_mask(const unsigned char *src, unsigned char *dst,
                        int width, int height) {
    memset(dst, 0, (size_t)width * height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (!raw_mask_at(src, width, height, x, y)) continue;
            for (int oy = -1; oy <= 1; ++oy) {
                int yy = y + oy;
                if (yy < 0 || yy >= height) continue;
                for (int ox = -1; ox <= 1; ++ox) {
                    int xx = x + ox;
                    if (xx < 0 || xx >= width) continue;
                    dst[(size_t)yy * width + xx] = 1;
                }
            }
        }
    }
}

static void erode_mask(const unsigned char *src, unsigned char *dst, int width,
                       int height) {
    memset(dst, 0, (size_t)width * height);
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int keep = 1;
            for (int oy = -1; oy <= 1 && keep; ++oy) {
                for (int ox = -1; ox <= 1; ++ox) {
                    if (!raw_mask_at(src, width, height, x + ox, y + oy)) {
                        keep = 0;
                        break;
                    }
                }
            }
            if (keep) dst[(size_t)y * width + x] = 1;
        }
    }
}

static void bridge_boundary_gaps(unsigned char *mask, int width, int height) {
    unsigned char *add = calloc((size_t)width * height, 1);
    if (add == NULL) die("Out of memory.");
    static const int dirs[8][2] = {{1, 0},  {0, 1},  {1, 1},  {1, -1},
                                   {-1, 0}, {0, -1}, {-1, -1}, {-1, 1}};
    for (int y = 2; y < height - 2; ++y) {
        for (int x = 2; x < width - 2; ++x) {
            if (raw_mask_at(mask, width, height, x, y)) continue;
            for (int i = 0; i < 8; ++i) {
                int dx = dirs[i][0];
                int dy = dirs[i][1];
                if (raw_mask_at(mask, width, height, x - dx, y - dy) &&
                    raw_mask_at(mask, width, height, x + dx, y + dy)) {
                    add[(size_t)y * width + x] = 1;
                    break;
                }
                if (raw_mask_at(mask, width, height, x - 2 * dx,
                                y - 2 * dy) &&
                    raw_mask_at(mask, width, height, x + dx, y + dy)) {
                    add[(size_t)y * width + x] = 1;
                    add[(size_t)(y - dy) * width + x - dx] = 1;
                    break;
                }
            }
        }
    }
    for (int i = 0; i < width * height; ++i) {
        if (add[i]) mask[i] = 1;
    }
    free(add);
}

static void prune_boundary_speckles(unsigned char *mask, int width,
                                    int height) {
    unsigned char *remove = calloc((size_t)width * height, 1);
    if (remove == NULL) die("Out of memory.");
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            if (!raw_mask_at(mask, width, height, x, y)) continue;
            if (raw_mask_neighbors(mask, width, height, x, y) <= 1) {
                remove[(size_t)y * width + x] = 1;
            }
        }
    }
    for (int i = 0; i < width * height; ++i) {
        if (remove[i]) mask[i] = 0;
    }
    free(remove);
}

static unsigned char *build_orientation_boundary_mask(
    const Orientation *orientation, int width, int height) {
    unsigned char *mask = calloc((size_t)width * height, 1);
    unsigned char *tmp = calloc((size_t)width * height, 1);
    if (mask == NULL || tmp == NULL) die("Out of memory.");
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            if (is_orientation_boundary(orientation, width, height, x, y)) {
                mask[(size_t)y * width + x] = 1;
            }
        }
    }
    dilate_mask(mask, tmp, width, height);
    erode_mask(tmp, mask, width, height);
    bridge_boundary_gaps(mask, width, height);
    thin_mask(mask, width, height);
    prune_boundary_speckles(mask, width, height);
    bridge_boundary_gaps(mask, width, height);
    free(tmp);
    return mask;
}

static int boundary_mask_sample(const unsigned char *mask, int width,
                                int height, int x, int y, int radius) {
    if (mask == NULL) return 0;
    for (int oy = -radius; oy <= radius; ++oy) {
        for (int ox = -radius; ox <= radius; ++ox) {
            if (raw_mask_at(mask, width, height, x + ox, y + oy)) return 1;
        }
    }
    return 0;
}

static void thin_mask(unsigned char *mask, int width, int height) {
    unsigned char *remove = calloc((size_t)width * height, 1);
    if (remove == NULL) die("Out of memory.");
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int pass = 0; pass < 2; ++pass) {
            memset(remove, 0, (size_t)width * height);
            for (int y = 1; y < height - 1; ++y) {
                for (int x = 1; x < width - 1; ++x) {
                    if (!mask_at(mask, width, x, y)) continue;
                    int p2 = mask_at(mask, width, x, y - 1);
                    int p3 = mask_at(mask, width, x + 1, y - 1);
                    int p4 = mask_at(mask, width, x + 1, y);
                    int p5 = mask_at(mask, width, x + 1, y + 1);
                    int p6 = mask_at(mask, width, x, y + 1);
                    int p7 = mask_at(mask, width, x - 1, y + 1);
                    int p8 = mask_at(mask, width, x - 1, y);
                    int p9 = mask_at(mask, width, x - 1, y - 1);
                    int n = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
                    if (n < 2 || n > 6) continue;
                    int transitions = (!p2 && p3) + (!p3 && p4) +
                                      (!p4 && p5) + (!p5 && p6) +
                                      (!p6 && p7) + (!p7 && p8) +
                                      (!p8 && p9) + (!p9 && p2);
                    if (transitions != 1) continue;
                    if (pass == 0) {
                        if (p2 && p4 && p6) continue;
                        if (p4 && p6 && p8) continue;
                    } else {
                        if (p2 && p4 && p8) continue;
                        if (p2 && p6 && p8) continue;
                    }
                    remove[(size_t)y * width + x] = 1;
                }
            }
            for (int y = 1; y < height - 1; ++y) {
                for (int x = 1; x < width - 1; ++x) {
                    size_t idx = (size_t)y * width + x;
                    if (remove[idx]) {
                        mask[idx] = 0;
                        changed = 1;
                    }
                }
            }
        }
    }
    free(remove);
}

static void UNUSED_FUNCTION render_boundaries(Canvas *canvas,
                                              const Field *field,
                                              const Orientation *orientation,
                                              const unsigned char *boundary_mask,
                                              int x0, int y0, int w, int h) {
    (void)orientation;
    rect(canvas, x0, y0, w, h, 252, 252, 252);
    int img_x = x0, img_y = y0, img_w = w, img_h = h;
    for (int i = 0; i <= 4; ++i) {
        int xx = img_x + i * img_w / 4;
        int yy = img_y + i * img_h / 4;
        line(canvas, xx, img_y, xx, img_y + img_h, 190, 190, 190);
        line(canvas, img_x, yy, img_x + img_w, yy, 190, 190, 190);
    }
    for (int y = 1; y < img_h - 1; ++y) {
        int sy = field->height - 1 -
                 (int)((double)y / (img_h - 1) * (field->height - 1) + 0.5);
        for (int x = 1; x < img_w - 1; ++x) {
            int sx =
                (int)((double)x / (img_w - 1) * (field->width - 1) + 0.5);
            if (boundary_mask_sample(boundary_mask, field->width,
                                     field->height, sx, sy, 1)) {
                pixel(canvas, img_x + x, img_y + y, 226, 54, 48);
            }
        }
    }
    line(canvas, img_x, img_y, img_x + img_w - 1, img_y, 120, 120, 120);
    line(canvas, img_x, img_y + img_h - 1, img_x + img_w - 1,
         img_y + img_h - 1, 120, 120, 120);
    line(canvas, img_x, img_y, img_x, img_y + img_h - 1, 120, 120, 120);
    line(canvas, img_x + img_w - 1, img_y, img_x + img_w - 1,
         img_y + img_h - 1, 120, 120, 120);

    text(canvas, img_x + img_w / 2 - 8, img_y + img_h + 14, "N", 2, 40, 40,
         50);
    text(canvas, img_x + img_w / 2 + 10, img_y + img_h + 29, "x", 1, 40, 40,
         50);
    text(canvas, img_x - 42, img_y + img_h / 2 - 8, "N", 2, 40, 40, 50);
    text(canvas, img_x - 24, img_y + img_h / 2 + 7, "y", 1, 40, 40, 50);
    line(canvas, img_x - 18, img_y + img_h + 12, img_x + 22,
         img_y + img_h + 12, 50, 50, 50);
    line(canvas, img_x - 18, img_y + img_h + 12, img_x - 18,
         img_y + img_h - 28, 50, 50, 50);
    line(canvas, img_x + 22, img_y + img_h + 12, img_x + 14,
         img_y + img_h + 8, 50, 50, 50);
    line(canvas, img_x + 22, img_y + img_h + 12, img_x + 14,
         img_y + img_h + 16, 50, 50, 50);
    line(canvas, img_x - 18, img_y + img_h - 28, img_x - 22,
         img_y + img_h - 20, 50, 50, 50);
    line(canvas, img_x - 18, img_y + img_h - 28, img_x - 14,
         img_y + img_h - 20, 50, 50, 50);
    text(canvas, img_x + 30, img_y + img_h + 10, "x", 1, 40, 40, 50);
    text(canvas, img_x - 34, img_y + img_h - 48, "y", 1, 40, 40, 50);
}

static void arrow(Canvas *canvas, int x0, int y0, int x1, int y1) {
    line(canvas, x0, y0, x1, y1, 45, 45, 50);
    line(canvas, x1, y1, x1 - 10, y1 - 5, 45, 45, 50);
    line(canvas, x1, y1, x1 - 10, y1 + 5, 45, 45, 50);
}

static void UNUSED_FUNCTION frame_label(const char *filename, char *label,
                                        size_t size) {
    const char *p = strstr(filename, "-t-");
    if (p == NULL) p = strstr(filename, "_t_");
    if (p == NULL) {
        snprintf(label, size, "t");
        return;
    }
    p += 3;
    char number[32];
    int n = 0;
    while (*p >= '0' && *p <= '9' && n < (int)sizeof(number) - 1) {
        number[n++] = *p++;
    }
    number[n] = '\0';
    if (n == 0) {
        snprintf(label, size, "t");
    } else {
        snprintf(label, size, "t=%s", number);
    }
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s --width W --height H --output figure.png "
            "[--style pfc|schematic] "
            "[--plot-variant plateau|scatter] "
            "[--panel-b-csv stats.csv] "
            "snap1.dat snap2.dat snap3.dat snap4.dat relaxed.dat\n",
            program);
}

int main(int argc, char **argv) {
    int width = 288;
    int height = 192;
    int plot_variant = 0;
    int style_schematic = 0;
    const char *output = "paper-figure.png";
    const char *panel_b_csv = NULL;
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
        } else if (strcmp(argv[i], "--style") == 0 && i + 1 < argc) {
            const char *name = argv[++i];
            if (strcmp(name, "pfc") == 0) {
                style_schematic = 0;
            } else if (strcmp(name, "schematic") == 0) {
                style_schematic = 1;
            } else {
                usage(argv[0]);
                return EXIT_FAILURE;
            }
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
        } else if (strcmp(argv[i], "--panel-b-csv") == 0 && i + 1 < argc) {
            panel_b_csv = argv[++i];
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
    Orientation orientation =
        compute_orientation(&fields[4], ORIENTATION_RADIUS);
    smooth_orientation(&orientation, fields[4].width, fields[4].height,
                       ORIENTATION_SMOOTH_ITERATIONS);
    unsigned char *boundary_mask = build_orientation_boundary_mask(
        &orientation, fields[4].width, fields[4].height);
    double gray_low, gray_high;
    gray_range(fields, 4, &gray_low, &gray_high);
    PlotData panel_b_data = {0};
    if (panel_b_csv != NULL) {
        panel_b_data = read_plot_csv(panel_b_csv);
    } else {
        static PlotPoint plateau_points[] = {
            {0.3, 2.85}, {0.7, 2.82}, {1.1, 2.86}, {1.5, 2.91},
            {1.9, 3.02}, {2.3, 3.16}, {2.7, 3.30}, {3.1, 3.50},
            {3.5, 3.75}, {3.9, 4.05}, {4.3, 4.40}, {4.8, 5.05},
            {5.6, 5.65}, {7.0, 6.90}, {8.5, 8.45}, {10.2, 9.80},
            {11.0, 11.0}, {12.5, 12.2}, {13.2, 14.0}};
        static PlotPoint scatter_points[] = {
            {0.5, 2.72}, {1.0, 3.08}, {1.7, 2.96}, {2.4, 3.46},
            {3.0, 3.25}, {3.8, 4.26}, {4.7, 4.62}, {5.4, 5.15},
            {6.4, 6.12}, {7.5, 6.88}, {8.4, 8.05}, {9.7, 8.66},
            {10.8, 10.7}, {12.0, 11.4}, {13.3, 13.1}};
        panel_b_data.points = plot_variant ? scatter_points : plateau_points;
        panel_b_data.count =
            plot_variant ? (int)(sizeof(scatter_points) / sizeof(scatter_points[0]))
                         : (int)(sizeof(plateau_points) / sizeof(plateau_points[0]));
    }

    Canvas canvas;
    canvas_init(&canvas, 1200, 600);

    int top_x = 70, top_y = 34, panel_size = 180;
    int panel_w = panel_size, panel_h = panel_size, gap = 16;
    int orient_x = top_x + 4 * (panel_w + gap);
    int timeline_y = top_y + panel_h + 12;
    int timeline_label_y = timeline_y + 7;
    int relaxed_label_y = top_y + panel_h + 8;
    int plot_x = 90, lower_y = 310, plot_w = 296, plot_h = 220;
    int sheet_x = 430, sheet_y = 318, sheet_w = 430, sheet_h = 240;
    int gb_x = 930, gb_y = 310, gb_w = 228, gb_h = 228;

    text(&canvas, 26, 8, "(a)", 2, 45, 45, 55);
    text(&canvas, plot_x - 64, lower_y - 26, "(b)", 2, 45, 45, 55);
    text(&canvas, sheet_x - 18, lower_y - 26, "(c)", 2, 45, 45, 55);
    text(&canvas, gb_x - 46, gb_y - 22, "(d)", 2, 45, 45, 55);

    Scene scene = {0};
    if (style_schematic) {
        scene = scene_create(width, height, 190);
        for (int i = 0; i < 4; ++i) {
            render_scene_growth(&canvas, &scene, top_x + i * (panel_w + gap),
                                top_y, panel_w, panel_h, i);
        }
        render_scene_orientation(&canvas, &scene, orient_x, top_y, panel_w,
                                 panel_h);
    } else {
        for (int i = 0; i < 4; ++i) {
            render_gray(&canvas, &fields[i], top_x + i * (panel_w + gap),
                        top_y, panel_w, panel_h, gray_low, gray_high, i);
        }
        render_orientation(&canvas, &fields[4], &orientation, boundary_mask,
                           orient_x, top_y, panel_w, panel_h, 1);
    }
    render_vertical_colorbar(&canvas, orient_x + panel_w + 28, top_y + 24, 18,
                             132);
    arrow(&canvas, top_x + 80, timeline_y,
          top_x + 3 * (panel_w + gap) + 120, timeline_y);
    text(&canvas, 286, timeline_label_y, "PFC simulation of CVD graphene growth", 1,
         45, 45, 50);
    text(&canvas, orient_x + 35, relaxed_label_y, "Relaxed nc-graphene", 1,
         45, 45, 50);

    render_plot(&canvas, plot_x, lower_y, plot_w, plot_h, &panel_b_data);
    if (style_schematic) {
        render_scene_3d(&canvas, &scene, sheet_x, sheet_y);
        render_scene_boundaries(&canvas, &scene, gb_x, gb_y, gb_w, gb_h);
    } else {
        render_3d(&canvas, &fields[4], &orientation, sheet_x, sheet_y,
                  sheet_w, sheet_h);
        render_boundaries(&canvas, &fields[4], &orientation, boundary_mask,
                          gb_x, gb_y, gb_w, gb_h);
    }

    field_write_rgb_png(canvas.rgb, canvas.width, canvas.height, output);

    canvas_free(&canvas);
    if (panel_b_csv != NULL) free_plot_data(&panel_b_data);
    if (style_schematic) {
        scene_free(&scene);
    }
    free(boundary_mask);
    free_orientation(&orientation);
    for (int i = 0; i < 5; ++i) free_field(&fields[i]);
    return EXIT_SUCCESS;
}
