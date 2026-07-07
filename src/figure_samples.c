#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884

typedef struct {
    double x;
    double y;
    double theta;
    double radius_bias;
    double nucleus_radius;
    double nucleus_contrast;
} Grain;

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static double rnd(void) {
    return (double)rand() / (double)RAND_MAX;
}

static double periodic_delta(double a, double b, double length) {
    double d = a - b;
    if (d > 0.5 * length) d -= length;
    if (d < -0.5 * length) d += length;
    return d;
}

static double one_mode(double x, double y, double lattice, double theta) {
    double qx = 2.0 * PI / lattice;
    double qy = qx / sqrt(3.0);
    double c = cos(theta);
    double s = sin(theta);
    double xr = c * x - s * y;
    double yr = s * x + c * y;
    return cos(qx * xr) * cos(qy * yr) + 0.5 * cos(2.0 * qy * yr);
}

static void write_field(const char *filename, const double *field, int width,
                        int height) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) die("Unable to write sample field.");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(file, "%.9e\n", field[(size_t)y * width + x]);
        }
    }
    fclose(file);
    printf("Wrote %s\n", filename);
}

static void make_grains(Grain *grains, int count, int width, int height) {
    int cols = (int)ceil(sqrt((double)count * width / height));
    int rows = (count + cols - 1) / cols;
    int index = 0;
    for (int row = 0; row < rows && index < count; ++row) {
        for (int col = 0; col < cols && index < count; ++col) {
            double cx = (col + 0.5) * width / cols;
            double cy = (row + 0.5) * height / rows;
            grains[index].x = cx + (rnd() - 0.5) * width / cols * 0.65;
            grains[index].y = cy + (rnd() - 0.5) * height / rows * 0.65;
            grains[index].theta = -PI / 12.0 + rnd() * PI / 6.0;
            grains[index].radius_bias = 0.82 + 0.38 * rnd();
            grains[index].nucleus_radius = 1.4 + 2.5 * rnd();
            grains[index].nucleus_contrast = 0.40 + 0.28 * rnd();
            ++index;
        }
    }
}

static int nearest_two_grains(const Grain *grains, int count, double x, double y,
                              int width, int height, double *best_dist2,
                              double *second_dist2) {
    int best = 0;
    double best_d2 = 1e100;
    double second_d2 = 1e100;
    for (int i = 0; i < count; ++i) {
        double dx = periodic_delta(x, grains[i].x, width);
        double dy = periodic_delta(y, grains[i].y, height);
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) {
            second_d2 = best_d2;
            best_d2 = d2;
            best = i;
        } else if (d2 < second_d2) {
            second_d2 = d2;
        }
    }
    *best_dist2 = best_d2;
    *second_dist2 = second_d2;
    return best;
}

static void generate_frame(double *field, int width, int height,
                           const Grain *grains, int grain_count, int frame,
                           int frame_count) {
    double dx = 0.7;
    double dy = 0.7;
    double lattice = 7.2552;
    double min_side = width < height ? width : height;
    double stage_radius[5] = {
        0.025 * min_side, 0.062 * min_side, 0.095 * min_side,
        0.145 * min_side, 1.20 * min_side};
    double radius = stage_radius[frame < 5 ? frame : 4];

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double d2;
            double second_d2;
            int g = nearest_two_grains(grains, grain_count, x, y, width, height,
                                       &d2, &second_d2);
            double rx = periodic_delta(x, grains[g].x, width) * dx;
            double ry = periodic_delta(y, grains[g].y, height) * dy;
            double dist = sqrt(d2);
            double second_dist = sqrt(second_d2);
            double local_radius = radius * grains[g].radius_bias;
            double edge = (dist - local_radius) / 1.8;
            double mask = frame == frame_count - 1 ? 1.0 :
                          1.0 / (1.0 + exp(edge));
            double wave = one_mode(rx, ry, lattice, grains[g].theta);
            double value = 0.94;
            if (frame == 0) {
                double nr = grains[g].nucleus_radius;
                double nuclei = exp(-d2 / (2.0 * nr * nr));
                double grain_noise = 0.035 * sin(0.19 * x + 0.31 * y + 7.0 * g);
                value = 0.95 - grains[g].nucleus_contrast * nuclei *
                              (0.72 + 0.28 * wave) + grain_noise * nuclei;
            } else {
                double interior = 0.54 - 0.20 * wave;
                double background = 0.94 + 0.010 * sin(0.07 * x + 0.03 * y);
                value = background * (1.0 - mask) + interior * mask;
                if (frame >= 2 && frame < frame_count - 1) {
                    double seam_wobble =
                        0.95 * sin(0.090 * x + 0.047 * y + 1.7 * g) +
                        0.55 * sin(0.033 * x - 0.111 * y);
                    double seam = fabs(second_dist - dist + seam_wobble);
                    double neighbor_radius = radius * grains[g].radius_bias;
                    double covered = mask > 0.55 && second_dist < neighbor_radius + 3.5;
                    if (covered && seam < 1.05) {
                        double t = seam / 1.05;
                        double boundary = 0.68 + 0.06 * sin(0.21 * x + 0.13 * y);
                        value = value * t + boundary * (1.0 - t);
                    }
                }
                if (frame == 3 && mask > 0.58 && mask < 0.86) {
                    double edge_noise = 0.015 * sin(0.23 * x - 0.17 * y + 3.0 * g);
                    value = 0.78 + edge_noise;
                }
                if (frame == frame_count - 1) {
                    double sx = periodic_delta(x, grains[g].x, width);
                    double sy = periodic_delta(y, grains[g].y, height);
                    double wobble = 0.018 * sin(0.050 * sx + 0.037 * sy);
                    value = 0.55 - 0.19 * wave + wobble;
                    double seam_wobble =
                        0.70 * sin(0.075 * x + 0.041 * y + 1.9 * g);
                    if (fabs(second_dist - dist + seam_wobble) < 1.25) {
                        value = 0.70 + 0.03 * sin(0.13 * x + 0.17 * y);
                    }
                }
            }
            field[(size_t)y * width + x] = value;
        }
    }
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s [--width 288] [--height 192] [--grains 36] "
            "[--prefix discrete]\n",
            program);
}

int main(int argc, char **argv) {
    int width = 288;
    int height = 192;
    int grains_n = 36;
    const char *prefix = "discrete";
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            width = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--height") == 0 && i + 1 < argc) {
            height = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--grains") == 0 && i + 1 < argc) {
            grains_n = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--prefix") == 0 && i + 1 < argc) {
            prefix = argv[++i];
        } else {
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    if (width <= 16 || height <= 16 || grains_n <= 0) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    srand(20260707u);
    Grain *grains = malloc((size_t)grains_n * sizeof(*grains));
    double *field = malloc((size_t)width * height * sizeof(*field));
    if (grains == NULL || field == NULL) die("Out of memory.");
    make_grains(grains, grains_n, width, height);

    for (int frame = 0; frame < 5; ++frame) {
        char filename[320];
        snprintf(filename, sizeof(filename), "%s-t-%d.dat", prefix, frame);
        generate_frame(field, width, height, grains, grains_n, frame, 5);
        write_field(filename, field, width, height);
    }

    free(field);
    free(grains);
    return EXIT_SUCCESS;
}
