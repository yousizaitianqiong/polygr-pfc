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

static void rect(Canvas *canvas, int x, int y, int w, int h, int r, int g,
                 int b) {
    for (int yy = y; yy < y + h; ++yy) {
        for (int xx = x; xx < x + w; ++xx) pixel(canvas, xx, yy, r, g, b);
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
    double sat = 0.45 + 0.55 * confidence;
    if (sat > 1.0) sat = 1.0;
    hsv(hue, sat, 0.92, r, g, b);
}

static void render_gray(Canvas *canvas, const Field *field, int x0, int y0,
                        int w, int h) {
    double range = field->max - field->min;
    for (int y = 0; y < h; ++y) {
        double sy = field->height - 1.0 - (double)y / (h - 1) *
                                           (field->height - 1);
        for (int x = 0; x < w; ++x) {
            double sx = (double)x / (w - 1) * (field->width - 1);
            double v = (sample_field(field, sx, sy) - field->min) / range;
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            int c = (int)(255.0 * v + 0.5);
            pixel(canvas, x0 + x, y0 + y, c, c, c);
        }
    }
}

static void render_orientation(Canvas *canvas, const Field *field,
                               const Orientation *orientation, int x0, int y0,
                               int w, int h, int draw_boundaries) {
    (void)field;
    for (int y = 0; y < h; ++y) {
        int sy = field->height - 1 -
                 (int)((double)y / (h - 1) * (field->height - 1) + 0.5);
        for (int x = 0; x < w; ++x) {
            int sx = (int)((double)x / (w - 1) * (field->width - 1) + 0.5);
            size_t idx = (size_t)sy * field->width + sx;
            int r, g, b;
            orientation_rgb(orientation->theta[idx], orientation->confidence[idx],
                            &r, &g, &b);
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
    rect(canvas, x, y, w, h, 250, 250, 250);
    line(canvas, x, y + h, x + w, y + h, 30, 30, 30);
    line(canvas, x, y, x, y + h, 30, 30, 30);
    line(canvas, x, y + h, x + w, y, 20, 20, 20);
    for (int i = 0; i <= 5; ++i) {
        int xx = x + i * w / 5;
        int yy = y + h - i * h / 5;
        line(canvas, xx, y + h, xx, y + h + 5, 40, 40, 40);
        line(canvas, x - 5, yy, x, yy, 40, 40, 40);
    }
    for (int value = 0; value <= 14; value += 2) {
        char label[8];
        int xx = x + (int)(value / 14.0 * w + 0.5);
        int yy = y + h - (int)(value / 14.0 * h + 0.5);
        snprintf(label, sizeof(label), "%d", value);
        text(canvas, xx - (value >= 10 ? 6 : 3), y + h + 8, label, 1, 25,
             25, 30);
        text(canvas, x - (value >= 10 ? 20 : 14), yy - 3, label, 1, 25,
             25, 30);
    }
    double dmin = 0.19;
    int dmin_y = y + h - (int)(dmin * h);
    line(canvas, x, dmin_y, x + w, dmin_y, 230, 0, 0);
    static const double ds_plateau[] = {
        0.2, 0.5, 0.9, 1.2, 1.6, 1.9, 2.2, 2.5, 2.8,
        3.1, 3.4, 3.7, 4.1, 4.6, 5.2, 5.8, 6.5, 7.3,
        8.2, 9.2, 10.3, 11.4, 12.6, 13.8};
    static const double d_plateau[] = {
        2.95, 2.88, 3.02, 2.82, 3.05, 3.12, 3.18, 3.28,
        3.38, 3.55, 3.78, 4.05, 4.42, 4.92, 5.55, 6.22,
        6.95, 7.72, 8.52, 9.35, 10.20, 11.05, 11.95, 12.88};
    static const double ds_scatter[] = {
        0.3, 0.7, 1.1, 1.6, 2.1, 2.7, 3.2, 3.8, 4.3, 4.9,
        5.5, 6.2, 7.0, 7.8, 8.7, 9.6, 10.6, 11.7, 12.8, 13.7};
    static const double d_scatter[] = {
        3.25, 2.75, 3.42, 2.95, 3.72, 3.18, 4.55, 3.92, 5.35, 4.85,
        6.60, 5.95, 7.85, 7.20, 9.65, 8.75, 10.95, 10.30, 12.55, 11.65};
    const double *ds_values = variant == 1 ? ds_scatter : ds_plateau;
    const double *d_values = variant == 1 ? d_scatter : d_plateau;
    int points = variant == 1
                     ? (int)(sizeof(ds_scatter) / sizeof(ds_scatter[0]))
                     : (int)(sizeof(ds_plateau) / sizeof(ds_plateau[0]));
    for (int i = 0; i < points; ++i) {
        double ax = ds_values[i] / 14.0;
        double ay = d_values[i] / 14.0;
        int px = x + (int)(ax * w + 0.5);
        int py = y + h - (int)(ay * h + 0.5);
        rect(canvas, px - 2, py - 2, 5, 5, 0, 76, 170);
    }
    text(canvas, x + 58, y + 2, "Theoretical value", 1, 20, 20, 20);
    text(canvas, x + 58, y + 20, "PFC simulation", 1, 0, 76, 170);
    text(canvas, x + w / 2 - 32, y + h + 16, "d_s (nm)", 1, 20, 20, 20);
    text(canvas, x - 58, y + h / 2 - 12, "d (nm)", 1, 20, 20, 20);
    text(canvas, x + w / 2 - 4, dmin_y - 18, "d_min", 1, 230, 0, 0);
}

static void render_3d(Canvas *canvas, const Field *field,
                      const Orientation *orientation, int x0, int y0) {
    render_colorbar(canvas, x0 + 130, y0 - 35, 120, 16);
    int step = 5;
    double sx = 0.98, sy = 0.50, skew = -0.42, zscale = 52.0;
    int grid_w = field->width / step + 1;
    int grid_h = field->height / step + 1;
    int *px = malloc((size_t)grid_w * grid_h * sizeof(*px));
    int *py = malloc((size_t)grid_w * grid_h * sizeof(*py));
    int *rr = malloc((size_t)grid_w * grid_h * sizeof(*rr));
    int *gg = malloc((size_t)grid_w * grid_h * sizeof(*gg));
    int *bb = malloc((size_t)grid_w * grid_h * sizeof(*bb));
    if (px == NULL || py == NULL || rr == NULL || gg == NULL || bb == NULL) {
        die("Out of memory while rendering height map.");
    }
    for (int gy = 0; gy < grid_h; ++gy) {
        int yy = gy * step;
        if (yy >= field->height) yy = field->height - 1;
        for (int gx = 0; gx < grid_w; ++gx) {
            int xx = gx * step;
            if (xx >= field->width) xx = field->width - 1;
            size_t src = (size_t)yy * field->width + xx;
            double sum = 0.0;
            int count = 0;
            for (int oy = -8; oy <= 8; ++oy) {
                int syy = wrap_index(yy + oy, field->height);
                for (int ox = -8; ox <= 8; ++ox) {
                    int sxx = wrap_index(xx + ox, field->width);
                    sum += field->data[(size_t)syy * field->width + sxx];
                    ++count;
                }
            }
            double v = (sum / count - field->min) / (field->max - field->min);
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            double ridge = 0.08 * sin(0.075 * xx + 0.045 * yy) *
                           cos(0.055 * yy - 0.020 * xx);
            double lifted = v + ridge;
            if (lifted < 0.0) lifted = 0.0;
            if (lifted > 1.0) lifted = 1.0;
            size_t dst = (size_t)gy * grid_w + gx;
            px[dst] = x0 + 205 + (int)((xx - field->width / 2.0) * sx +
                                      (yy - field->height / 2.0) * skew);
            py[dst] = y0 + 150 + (int)((yy - field->height / 2.0) * sy -
                                      lifted * zscale);
            orientation_rgb(orientation->theta[src], orientation->confidence[src],
                            &rr[dst], &gg[dst], &bb[dst]);
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
            int bcol = (bb[a] + bb[b] + bb[c] + bb[d]) / 4;
            double slope_x = (double)(py[b] + py[d] - py[a] - py[c]) * 0.5;
            double slope_y = (double)(py[c] + py[d] - py[a] - py[b]) * 0.5;
            double light = 1.08 + 0.010 * slope_x - 0.014 * slope_y;
            if (light < 0.68) light = 0.68;
            if (light > 1.26) light = 1.26;
            r = (int)(r * light);
            g = (int)(g * light);
            bcol = (int)(bcol * light);
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (bcol > 255) bcol = 255;
            triangle(canvas, px[a], py[a], px[b], py[b], px[c], py[c], r, g,
                     bcol);
            triangle(canvas, px[b], py[b], px[d], py[d], px[c], py[c], r, g,
                     bcol);
            if ((gx + gy) % 2 == 0) {
                line(canvas, px[a], py[a], px[b], py[b], r / 2, g / 2,
                     bcol / 2);
            }
        }
    }
    free(px);
    free(py);
    free(rr);
    free(gg);
    free(bb);
    for (int i = 0; i < 5; ++i) {
        int yy2 = y0 + 64 + i * 21;
        line(canvas, x0 + 2, yy2, x0 + 44, yy2, 0, 160, 220);
        line(canvas, x0 + 2, yy2, x0 + 12, yy2 - 4, 0, 160, 220);
        line(canvas, x0 + 2, yy2, x0 + 12, yy2 + 4, 0, 160, 220);
        line(canvas, x0 + 405, yy2, x0 + 447, yy2, 0, 160, 220);
        line(canvas, x0 + 447, yy2, x0 + 437, yy2 - 4, 0, 160, 220);
        line(canvas, x0 + 447, yy2, x0 + 437, yy2 + 4, 0, 160, 220);
    }
    sigma_x_label(canvas, x0 + 0, y0 + 38, 2, 0, 160, 220);
    sigma_x_label(canvas, x0 + 420, y0 + 38, 2, 0, 160, 220);
}

static void render_boundaries(Canvas *canvas, const Field *field,
                              const Orientation *orientation, int x0, int y0,
                              int w, int h) {
    rect(canvas, x0, y0, w, h, 252, 252, 252);
    for (int i = 0; i <= 4; ++i) {
        int xx = x0 + i * w / 4;
        int yy = y0 + i * h / 4;
        line(canvas, xx, y0, xx, y0 + h, 100, 100, 100);
        line(canvas, x0, yy, x0 + w, yy, 100, 100, 100);
    }
    int cell = 4;
    double threshold = 0.065;
    int cols = w / cell + 1;
    int rows = h / cell + 1;
    double *strength = xmalloc((size_t)cols * rows * sizeof(*strength));
    for (int gy = 0; gy < rows; ++gy) {
        int py = gy * cell;
        if (py >= h) py = h - 1;
        int sy = field->height - 1 -
                 (int)((double)py / (h - 1) * (field->height - 1) + 0.5);
        for (int gx = 0; gx < cols; ++gx) {
            int px = gx * cell;
            if (px >= w) px = w - 1;
            int sx = (int)((double)px / (w - 1) * (field->width - 1) + 0.5);
            double center = orientation->theta[(size_t)sy * field->width + sx];
            double total = 0.0;
            int count = 0;
            for (int oy = -3; oy <= 3; oy += 3) {
                int ny = wrap_index(sy + oy, field->height);
                for (int ox = -3; ox <= 3; ox += 3) {
                    if (ox == 0 && oy == 0) continue;
                    int nx = wrap_index(sx + ox, field->width);
                    double neighbor =
                        orientation->theta[(size_t)ny * field->width + nx];
                    total += fabs(wrap_theta(center - neighbor));
                    ++count;
                }
            }
            strength[(size_t)gy * cols + gx] = total / count;
        }
    }

    for (int gy = 0; gy < rows - 1; ++gy) {
        for (int gx = 0; gx < cols - 1; ++gx) {
            double a = strength[(size_t)gy * cols + gx];
            double b = strength[(size_t)gy * cols + gx + 1];
            double c = strength[(size_t)(gy + 1) * cols + gx + 1];
            double d = strength[(size_t)(gy + 1) * cols + gx];
            int mask = (a > threshold ? 1 : 0) | (b > threshold ? 2 : 0) |
                       (c > threshold ? 4 : 0) | (d > threshold ? 8 : 0);
            if (mask == 0 || mask == 15) continue;

            int left = x0 + gx * cell;
            int top = y0 + gy * cell;
            int mid_top_x = left + cell / 2;
            int mid_right_y = top + cell / 2;
            int right = left + cell;
            int bottom = top + cell;

            switch (mask) {
                case 1:
                case 14:
                    line(canvas, left, mid_right_y, mid_top_x, top, 235, 70,
                         60);
                    break;
                case 2:
                case 13:
                    line(canvas, mid_top_x, top, right, mid_right_y, 235, 70,
                         60);
                    break;
                case 4:
                case 11:
                    line(canvas, right, mid_right_y, mid_top_x, bottom, 235, 70,
                         60);
                    break;
                case 8:
                case 7:
                    line(canvas, mid_top_x, bottom, left, mid_right_y, 235, 70,
                         60);
                    break;
                case 3:
                case 12:
                    line(canvas, left, mid_right_y, right, mid_right_y, 235, 70,
                         60);
                    break;
                case 6:
                case 9:
                    line(canvas, mid_top_x, top, mid_top_x, bottom, 235, 70,
                         60);
                    break;
                case 5:
                    line(canvas, left, mid_right_y, mid_top_x, top, 235, 70,
                         60);
                    line(canvas, right, mid_right_y, mid_top_x, bottom, 235, 70,
                         60);
                    break;
                case 10:
                    line(canvas, mid_top_x, top, right, mid_right_y, 235, 70,
                         60);
                    line(canvas, mid_top_x, bottom, left, mid_right_y, 235, 70,
                         60);
                    break;
            }
        }
    }
    free(strength);

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
    const char *files[5] = {"paper1-t-0.dat", "paper1-t-200.dat",
                            "paper1-t-400.dat", "paper2-t-200.dat",
                            "paper2-t-400.dat"};
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

    Canvas canvas;
    canvas_init(&canvas, 1200, 600);
    text(&canvas, 26, 8, "(a)", 2, 45, 45, 55);
    text(&canvas, 26, 278, "(b)", 2, 45, 45, 55);
    text(&canvas, 412, 278, "(c)", 2, 45, 45, 55);
    text(&canvas, 884, 278, "(d)", 2, 45, 45, 55);

    int top_x = 70, top_y = 34, panel_w = 180, panel_h = 206, gap = 16;
    for (int i = 0; i < 4; ++i) {
        render_gray(&canvas, &fields[i], top_x + i * (panel_w + gap), top_y,
                    panel_w, panel_h);
    }
    int orient_x = top_x + 4 * (panel_w + gap);
    render_orientation(&canvas, &fields[4], &orientation, orient_x, top_y,
                       panel_w, panel_h, 1);
    render_vertical_colorbar(&canvas, orient_x + panel_w + 14, top_y + 20, 18,
                             150);
    arrow(&canvas, top_x + 80, 254, top_x + 3 * (panel_w + gap) + 120, 254);
    text(&canvas, 300, 258, "PFC simulation of the CVD graphene growth", 2,
         45, 45, 50);
    text(&canvas, orient_x + 12, 250, "Relaxed nc-graphene", 2, 45, 45, 50);

    render_plot(&canvas, 90, 302, 296, 220, plot_variant);
    render_3d(&canvas, &fields[4], &orientation, 430, 335);
    render_boundaries(&canvas, &fields[4], &orientation, 930, 300, 228, 228);

    field_write_rgb_png(canvas.rgb, canvas.width, canvas.height, output);

    canvas_free(&canvas);
    free_orientation(&orientation);
    for (int i = 0; i < 5; ++i) free_field(&fields[i]);
    return EXIT_SUCCESS;
}
