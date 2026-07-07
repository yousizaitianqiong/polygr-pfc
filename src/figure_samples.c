#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884

typedef struct {
    double x;
    double y;
    double theta;
    double radius;
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
            grains[index].radius = 8.0 + 10.0 * rnd();
            ++index;
        }
    }
}

static int nearest_grain(const Grain *grains, int count, double x, double y,
                         int width, int height, double *dist2) {
    int best = 0;
    double best_d2 = 1e100;
    for (int i = 0; i < count; ++i) {
        double dx = periodic_delta(x, grains[i].x, width);
        double dy = periodic_delta(y, grains[i].y, height);
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) {
            best_d2 = d2;
            best = i;
        }
    }
    *dist2 = best_d2;
    return best;
}

static void generate_frame(double *field, int width, int height,
                           const Grain *grains, int grain_count, int frame,
                           int frame_count) {
    double dx = 0.7;
    double dy = 0.7;
    double lattice = 7.2552;
    double progress = (double)frame / (double)(frame_count - 1);
    double growth = 5.0 + progress * progress * 58.0;
    double density = 0.22;
    double amplitude = -0.18;

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double d2;
            int g = nearest_grain(grains, grain_count, x, y, width, height, &d2);
            double rx = periodic_delta(x, grains[g].x, width) * dx;
            double ry = periodic_delta(y, grains[g].y, height) * dy;
            double mask;
            if (frame == frame_count - 1) {
                mask = 1.0;
            } else {
                double local_radius = grains[g].radius + growth;
                double edge = (sqrt(d2) - local_radius) / 2.0;
                mask = 1.0 / (1.0 + exp(edge));
            }

            double wave = one_mode(rx, ry, lattice, grains[g].theta);
            double value = density + amplitude * mask * wave;
            if (frame == 0) {
                double nuclei = exp(-d2 / 14.0);
                value = 0.34 - 0.32 * nuclei * (0.55 + 0.45 * wave);
            }
            if (mask < 0.02 && frame > 0) {
                value = 0.34 + 0.004 * sin(0.37 * x + 0.19 * y);
            }
            field[(size_t)y * width + x] = value;
        }
    }
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s [--width 288] [--height 192] [--grains 72] "
            "[--prefix discrete]\n",
            program);
}

int main(int argc, char **argv) {
    int width = 288;
    int height = 192;
    int grains_n = 72;
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
