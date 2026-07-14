#include "field_image.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884

typedef struct {
    int width;
    int height;
    double *data;
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

static double wrap_theta(double theta) {
    double period = PI / 3.0;
    theta = fmod(theta + PI / 6.0, period);
    if (theta < 0.0) theta += period;
    return theta - PI / 6.0;
}

static Field read_field(const char *filename, int width, int height) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) die("Unable to open density field.");
    Field field = {width, height, xmalloc((size_t)width * height *
                                          sizeof(*field.data))};
    for (size_t i = 0; i < (size_t)width * height; ++i) {
        if (fscanf(file, "%lf", &field.data[i]) != 1) {
            die("Invalid density field data.");
        }
    }
    fclose(file);
    return field;
}

static void free_field(Field *field) {
    free(field->data);
    memset(field, 0, sizeof(*field));
}

static Orientation compute_orientation(const Field *field, int radius) {
    Orientation out = {xmalloc((size_t)field->width * field->height *
                               sizeof(*out.theta)),
                       xmalloc((size_t)field->width * field->height *
                               sizeof(*out.confidence))};
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
    double *btheta = xmalloc((size_t)bw * bh * sizeof(*btheta));
    double *bconf = xmalloc((size_t)bw * bh * sizeof(*bconf));

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
    double *next_conf = xmalloc((size_t)width * height * sizeof(*next_conf));
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
                        if (fabs(wrap_theta(theta - center)) > 0.070) {
                            continue;
                        }
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

static void pixel(unsigned char *rgb, int width, int height, int x, int y,
                  int r, int g, int b) {
    if (x < 0 || y < 0 || x >= width || y >= height) return;
    unsigned char *p = rgb + ((size_t)y * width + x) * 3;
    p[0] = (unsigned char)r;
    p[1] = (unsigned char)g;
    p[2] = (unsigned char)b;
}

static void line(unsigned char *rgb, int width, int height, int x0, int y0,
                 int x1, int y1, int r, int g, int b) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        pixel(rgb, width, height, x0, y0, r, g, b);
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

static double local_theta(const Orientation *orientation, int width, int height,
                          int x, int y, int radius) {
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

static void draw_grid(unsigned char *rgb, int width, int height, int margin,
                      int nx, int ny) {
    for (int i = 0; i <= nx; ++i) {
        int x = margin + (width - 2 * margin) * i / nx;
        line(rgb, width, height, x, margin, x, height - margin, 205, 205, 205);
    }
    for (int i = 0; i <= ny; ++i) {
        int y = margin + (height - 2 * margin) * i / ny;
        line(rgb, width, height, margin, y, width - margin, y, 205, 205, 205);
    }
}

static int mask_at(const unsigned char *mask, int width, int x, int y) {
    return mask[(size_t)y * width + x] != 0;
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

static void write_orientation_gbmap(const Field *field,
                                    const Orientation *orientation,
                                    const char *output, int image_w,
                                    int image_h, int grid_x, int grid_y,
                                    double threshold) {
    int margin = 28;
    unsigned char *rgb = xmalloc((size_t)image_w * image_h * 3);
    unsigned char *mask = calloc((size_t)image_w * image_h, 1);
    if (mask == NULL) die("Out of memory.");
    memset(rgb, 255, (size_t)image_w * image_h * 3);
    if (grid_x > 0 && grid_y > 0) {
        draw_grid(rgb, image_w, image_h, margin, grid_x, grid_y);
    }

    int count = 0;
    for (int y = margin; y < image_h - margin; ++y) {
        int sy = field->height - 1 -
                 (int)((double)(y - margin) /
                           (image_h - 2 * margin - 1) *
                           (field->height - 1) +
                       0.5);
        for (int x = margin; x < image_w - margin; ++x) {
            int sx = (int)((double)(x - margin) /
                               (image_w - 2 * margin - 1) *
                               (field->width - 1) +
                           0.5);
            double center = local_theta(orientation, field->width,
                                        field->height, sx, sy, 4);
            double score = 0.0;
            const int offsets[4][2] = {{6, 0}, {-6, 0}, {0, 6}, {0, -6}};
            for (int i = 0; i < 4; ++i) {
                int xx = wrap_index(sx + offsets[i][0], field->width);
                int yy = wrap_index(sy + offsets[i][1], field->height);
                double neighbor = local_theta(orientation, field->width,
                                              field->height, xx, yy, 4);
                double delta = fabs(wrap_theta(center - neighbor));
                if (delta > score) score = delta;
            }
            if (score <= threshold) continue;
            int neighbors = 0;
            for (int oy = -3; oy <= 3; ++oy) {
                int yy = wrap_index(sy + oy, field->height);
                for (int ox = -3; ox <= 3; ++ox) {
                    int xx = wrap_index(sx + ox, field->width);
                    double neighbor = local_theta(orientation, field->width,
                                                  field->height, xx, yy, 4);
                    if (fabs(wrap_theta(center - neighbor)) >
                        threshold * 0.78) {
                        ++neighbors;
                    }
                }
            }
            if (neighbors < 5) continue;
            mask[(size_t)y * image_w + x] = 1;
            ++count;
        }
    }

    thin_mask(mask, image_w, image_h);

    for (int y = margin; y < image_h - margin; ++y) {
        for (int x = margin; x < image_w - margin; ++x) {
            if (!mask[(size_t)y * image_w + x]) continue;
            pixel(rgb, image_w, image_h, x, y, 238, 72, 61);
        }
    }

    field_write_rgb_png(rgb, image_w, image_h, output);
    printf("Wrote orientation grain-boundary map with %d boundary pixels to %s\n",
           count, output);
    free(mask);
    free(rgb);
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s field.dat --width W --height H --output gb.png "
            "[--image-width 768] [--image-height 512] [--grid NX NY] "
            "[--threshold 0.135] [--radius 9]\n",
            program);
}

int main(int argc, char **argv) {
    const char *input = NULL;
    const char *output = "orientation-gb-map.png";
    int width = 0;
    int height = 0;
    int image_w = 768;
    int image_h = 512;
    int grid_x = 4;
    int grid_y = 4;
    int radius = 9;
    double threshold = 0.135;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            width = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--height") == 0 && i + 1 < argc) {
            height = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output = argv[++i];
        } else if (strcmp(argv[i], "--image-width") == 0 && i + 1 < argc) {
            image_w = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--image-height") == 0 && i + 1 < argc) {
            image_h = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--grid") == 0 && i + 2 < argc) {
            grid_x = atoi(argv[++i]);
            grid_y = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--threshold") == 0 && i + 1 < argc) {
            threshold = atof(argv[++i]);
        } else if (strcmp(argv[i], "--radius") == 0 && i + 1 < argc) {
            radius = atoi(argv[++i]);
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            return EXIT_FAILURE;
        } else if (input == NULL) {
            input = argv[i];
        } else {
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (input == NULL || width < 2 || height < 2 || image_w < 64 ||
        image_h < 64 || grid_x < 0 || grid_y < 0 || radius < 2 ||
        threshold <= 0.0) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    Field field = read_field(input, width, height);
    Orientation orientation = compute_orientation(&field, radius);
    smooth_orientation(&orientation, width, height, 2);
    write_orientation_gbmap(&field, &orientation, output, image_w, image_h,
                            grid_x, grid_y, threshold);
    free_orientation(&orientation);
    free_field(&field);
    return EXIT_SUCCESS;
}
