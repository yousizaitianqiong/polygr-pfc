/*
 * Pure C, shared-memory implementation of the polycrystalline graphene PFC
 * workflow. The phase-field solver uses FFTW threads and OpenMP. The final
 * density field is converted directly to XYZ without the original Java tools.
 */

#include <fftw3.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PI 3.141592653589793238462643383279502884

typedef struct {
    char name[128];
    int print_interval;
    int write_interval;
} Output;

typedef struct {
    int w;
    int h;
    int wc;
    size_t real_count;
    size_t complex_count;
    double *linear;
    double *nonlinear;
    double *p;
    double *q;
    fftw_plan p_forward;
    fftw_plan q_forward;
    fftw_plan q_inverse;
} Arrays;

typedef struct {
    double alpha;
    double beta;
    double gamma;
    double delta;
    double conserved;
} Model;

typedef struct {
    time_t started;
    int step;
    int total_steps;
    int optimize_interval;
    double free_energy;
    double density;
    double sample_step;
    double dx;
    double dy;
    double dt;
} Relaxation;

typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    int *ids;
    int count;
    int capacity;
} Bin;

typedef struct {
    Point *items;
    size_t count;
    size_t capacity;
} PointVector;

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void *checked_malloc(size_t bytes) {
    void *ptr = malloc(bytes);
    if (ptr == NULL) {
        die("Out of memory.");
    }
    return ptr;
}

static void point_vector_push(PointVector *vector, Point point) {
    if (vector->count == vector->capacity) {
        size_t next = vector->capacity == 0 ? 1024 : vector->capacity * 2;
        Point *items = realloc(vector->items, next * sizeof(*items));
        if (items == NULL) {
            die("Out of memory while storing carbon coordinates.");
        }
        vector->items = items;
        vector->capacity = next;
    }
    vector->items[vector->count++] = point;
}

static void bin_push(Bin *bin, int id) {
    if (bin->count == bin->capacity) {
        int next = bin->capacity == 0 ? 8 : bin->capacity * 2;
        int *ids = realloc(bin->ids, (size_t)next * sizeof(*ids));
        if (ids == NULL) {
            die("Out of memory while building spatial bins.");
        }
        bin->ids = ids;
        bin->capacity = next;
    }
    bin->ids[bin->count++] = id;
}

static int wrap_index(int value, int length) {
    value %= length;
    return value < 0 ? value + length : value;
}

static double wrap_coordinate(double value, double length) {
    value = fmod(value, length);
    return value < 0.0 ? value + length : value;
}

static Point periodic_difference(Point b, Point a, double width, double height) {
    Point result = {b.x - a.x, b.y - a.y};
    if (result.x > 0.5 * width) {
        result.x -= width;
    } else if (result.x < -0.5 * width) {
        result.x += width;
    }
    if (result.y > 0.5 * height) {
        result.y -= height;
    } else if (result.y < -0.5 * height) {
        result.y += height;
    }
    return result;
}

static double point_norm2(Point p) {
    return p.x * p.x + p.y * p.y;
}

static void configure_arrays(Arrays *arrays, FILE *input) {
    if (fscanf(input, " %d %d", &arrays->w, &arrays->h) != 2 ||
        arrays->w < 2 || arrays->h < 2) {
        die("Invalid array dimensions.");
    }

    arrays->wc = arrays->w / 2 + 1;
    arrays->real_count = (size_t)arrays->h * (size_t)(2 * arrays->wc);
    arrays->complex_count = (size_t)arrays->h * (size_t)arrays->wc;
    arrays->linear = fftw_alloc_real(arrays->complex_count);
    arrays->nonlinear = fftw_alloc_real(arrays->complex_count);
    arrays->p = fftw_alloc_real(arrays->real_count);
    arrays->q = fftw_alloc_real(arrays->real_count);
    if (!arrays->linear || !arrays->nonlinear || !arrays->p || !arrays->q) {
        die("FFTW allocation failed.");
    }

    arrays->p_forward = fftw_plan_dft_r2c_2d(
        arrays->h, arrays->w, arrays->p, (fftw_complex *)arrays->p, FFTW_MEASURE);
    arrays->q_forward = fftw_plan_dft_r2c_2d(
        arrays->h, arrays->w, arrays->q, (fftw_complex *)arrays->q, FFTW_MEASURE);
    arrays->q_inverse = fftw_plan_dft_c2r_2d(
        arrays->h, arrays->w, (fftw_complex *)arrays->q, arrays->q, FFTW_MEASURE);
    if (!arrays->p_forward || !arrays->q_forward || !arrays->q_inverse) {
        die("FFTW plan creation failed.");
    }
}

static void clear_arrays(Arrays *arrays) {
    fftw_destroy_plan(arrays->p_forward);
    fftw_destroy_plan(arrays->q_forward);
    fftw_destroy_plan(arrays->q_inverse);
    fftw_free(arrays->p);
    fftw_free(arrays->q);
    fftw_free(arrays->linear);
    fftw_free(arrays->nonlinear);
    memset(arrays, 0, sizeof(*arrays));
}

static void seed_rng(FILE *input) {
    unsigned int seed;
    if (fscanf(input, " %u", &seed) != 1) {
        die("Invalid random seed.");
    }
    srand(seed == 0 ? (unsigned int)time(NULL) : seed);
}

static void configure_output(Output *output, FILE *input) {
    if (fscanf(input, " %d %d", &output->print_interval,
               &output->write_interval) != 2 ||
        output->print_interval <= 0 || output->write_interval <= 0) {
        die("Invalid output configuration.");
    }
}

static void configure_model(Model *model, FILE *input) {
    if (fscanf(input, " %lf %lf %lf %lf %lf", &model->alpha, &model->beta,
               &model->gamma, &model->delta, &model->conserved) != 5) {
        die("Invalid model parameters.");
    }
}

static void configure_relaxation(Relaxation *relaxation, FILE *input) {
    if (fscanf(input, " %d %lf %lf %lf %d", &relaxation->total_steps,
               &relaxation->dx, &relaxation->dy, &relaxation->dt,
               &relaxation->optimize_interval) != 5) {
        die("Invalid relaxation settings.");
    }
}

static double one_mode(double x0, double y0, double lattice, double theta) {
    double qx = 2.0 * PI / lattice;
    double qy = qx / sqrt(3.0);
    double cosine = cos(theta);
    double sine = sin(theta);
    double x = cosine * x0 - sine * y0;
    double y = sine * x0 + cosine * y0;
    return cos(qx * x) * cos(qy * y) + 0.5 * cos(2.0 * qy * y);
}

static void random_state(Arrays *arrays, double density, double amplitude) {
    int stride = 2 * arrays->wc;
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            arrays->q[(size_t)y * stride + x] =
                density + amplitude * ((double)rand() / RAND_MAX - 0.5);
        }
    }
}

static void polycrystalline_state(Arrays *arrays, double dx, double dy,
                                  double lattice, double density,
                                  double amplitude, int grain_count,
                                  double radius) {
    if (grain_count <= 0 || radius <= 0.0) {
        die("Invalid polycrystalline initialization.");
    }
    double *grains =
        checked_malloc((size_t)grain_count * 3 * sizeof(*grains));
    for (int grain = 0; grain < grain_count; ++grain) {
        grains[3 * grain] = (double)rand() / RAND_MAX * arrays->w;
        grains[3 * grain + 1] = (double)rand() / RAND_MAX * arrays->h;
        grains[3 * grain + 2] = 2.0 * PI * rand() / RAND_MAX;
    }

    int stride = 2 * arrays->wc;
    double divisor = 0.5 / (radius * radius * lattice * lattice);
#pragma omp parallel for schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            int closest = 0;
            double closest_r2 = HUGE_VAL;
            double closest_x = 0.0;
            double closest_y = 0.0;
            for (int grain = 0; grain < grain_count; ++grain) {
                double rx = x - grains[3 * grain];
                double ry = y - grains[3 * grain + 1];
                if (rx > 0.5 * arrays->w) rx -= arrays->w;
                if (rx < -0.5 * arrays->w) rx += arrays->w;
                if (ry > 0.5 * arrays->h) ry -= arrays->h;
                if (ry < -0.5 * arrays->h) ry += arrays->h;
                rx *= dx;
                ry *= dy;
                double r2 = rx * rx + ry * ry;
                if (r2 < closest_r2) {
                    closest = grain;
                    closest_r2 = r2;
                    closest_x = rx;
                    closest_y = ry;
                }
            }
            arrays->q[(size_t)y * stride + x] =
                density + amplitude * exp(-closest_r2 * divisor) *
                              one_mode(closest_x, closest_y, lattice,
                                       grains[3 * closest + 2]);
        }
    }
    free(grains);
}

static void read_state(Arrays *arrays, FILE *file, double density,
                       double amplitude) {
    int stride = 2 * arrays->wc;
    double mean = 0.0;
    double minimum = HUGE_VAL;
    double maximum = -HUGE_VAL;
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            double value;
            if (fscanf(file, "%lf", &value) != 1) {
                die("Invalid state data.");
            }
            arrays->q[(size_t)y * stride + x] = value;
            mean += value;
            if (value < minimum) minimum = value;
            if (value > maximum) maximum = value;
        }
    }
    mean /= (double)arrays->w * arrays->h;
    double range = maximum - minimum;
    if (range == 0.0) {
        die("Cannot normalize a constant input state.");
    }
#pragma omp parallel for schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            size_t index = (size_t)y * stride + x;
            arrays->q[index] =
                (arrays->q[index] - mean) / range * amplitude + density;
        }
    }
}

static void initialize_system(Arrays *arrays, FILE *input) {
    int type;
    if (fscanf(input, " %d", &type) != 1) {
        die("Invalid initialization type.");
    }
    if (type == 0) {
        double density, amplitude;
        if (fscanf(input, " %lf %lf", &density, &amplitude) != 2) {
            die("Invalid random initialization.");
        }
        random_state(arrays, density, amplitude);
    } else if (type == 1) {
        int grains;
        double dx, dy, lattice, density, amplitude, radius;
        if (fscanf(input, " %lf %lf %lf %lf %lf %d %lf", &dx, &dy,
                   &lattice, &density, &amplitude, &grains, &radius) != 7) {
            die("Invalid polycrystalline initialization.");
        }
        polycrystalline_state(arrays, dx, dy, lattice, density, amplitude,
                              grains, radius);
    } else if (type == 2) {
        char filename[256];
        double density, amplitude;
        if (fscanf(input, " %255s %lf %lf", filename, &density, &amplitude) !=
            3) {
            die("Invalid file initialization.");
        }
        FILE *state = fopen(filename, "r");
        if (state == NULL) {
            fprintf(stderr, "State file not found: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        read_state(arrays, state, density, amplitude);
        fclose(state);
    } else {
        die("Unknown initialization type.");
    }
}

static void update_operators(Arrays *arrays, const Model *model,
                             const Relaxation *relaxation) {
    double inv_size = 1.0 / ((double)arrays->w * arrays->h);
    double dkx = 2.0 * PI / (relaxation->dx * arrays->w);
    double dky = 2.0 * PI / (relaxation->dy * arrays->h);
#pragma omp parallel for collapse(2) schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->wc; ++x) {
            double kx = x * dkx;
            double ky = (y < arrays->h / 2 ? y : y - arrays->h) * dky;
            double k2 = kx * kx + ky * ky;
            double one_minus_k2 = 1.0 - k2;
            double linear =
                model->alpha + model->beta * one_minus_k2 * one_minus_k2;
            double mobility = (1.0 - model->conserved) +
                              model->conserved * k2;
            double exponential = exp(-mobility * linear * relaxation->dt);
            size_t index = (size_t)y * arrays->wc + x;
            arrays->linear[index] = exponential;
            arrays->nonlinear[index] =
                (linear == 0.0 ? -mobility * relaxation->dt
                               : (exponential - 1.0) / linear) *
                inv_size;
        }
    }
}

static void scale_complex(double *data, size_t count, double scale) {
    fftw_complex *complex_data = (fftw_complex *)data;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count; ++i) {
        complex_data[i][0] *= scale;
        complex_data[i][1] *= scale;
    }
}

static void calculate_properties(Arrays *arrays, const Model *model,
                                 Relaxation *relaxation) {
    int stride = 2 * arrays->wc;
    memcpy(arrays->p, arrays->q, arrays->real_count * sizeof(*arrays->q));
    fftw_execute(arrays->q_forward);
    scale_complex(arrays->q, arrays->complex_count,
                  1.0 / ((double)arrays->w * arrays->h));

    double dkx = 2.0 * PI / (relaxation->dx * arrays->w);
    double dky = 2.0 * PI / (relaxation->dy * arrays->h);
    fftw_complex *q_frequency = (fftw_complex *)arrays->q;
#pragma omp parallel for collapse(2) schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->wc; ++x) {
            double kx = x * dkx;
            double ky = (y < arrays->h / 2 ? y : y - arrays->h) * dky;
            double factor = 1.0 - kx * kx - ky * ky;
            factor *= factor;
            size_t index = (size_t)y * arrays->wc + x;
            q_frequency[index][0] *= factor;
            q_frequency[index][1] *= factor;
        }
    }
    fftw_execute(arrays->q_inverse);

    double free_energy = 0.0;
    double density = 0.0;
#pragma omp parallel for reduction(+ : free_energy, density) schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            size_t index = (size_t)y * stride + x;
            double p = arrays->p[index];
            double p2 = p * p;
            free_energy +=
                0.5 * (model->alpha * p2 +
                       model->beta * p * arrays->q[index]) +
                model->gamma * p * p2 / 3.0 +
                0.25 * model->delta * p2 * p2;
            density += p;
            arrays->q[index] = p;
        }
    }
    double inv_size = 1.0 / ((double)arrays->w * arrays->h);
    relaxation->free_energy = free_energy * inv_size;
    relaxation->density = density * inv_size;
    fftw_execute(arrays->p_forward);
    scale_complex(arrays->p, arrays->complex_count, inv_size);
}

static void solver_step(Arrays *arrays, const Model *model,
                        const Relaxation *relaxation) {
    int stride = 2 * arrays->wc;
    if (relaxation->step % 100 == 0) {
        memcpy(arrays->p, arrays->q, arrays->real_count * sizeof(*arrays->q));
        fftw_execute(arrays->p_forward);
        scale_complex(arrays->p, arrays->complex_count,
                      1.0 / ((double)arrays->w * arrays->h));
    }

#pragma omp parallel for schedule(static)
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            size_t index = (size_t)y * stride + x;
            double value = arrays->q[index];
            arrays->q[index] =
                model->gamma * value * value +
                model->delta * value * value * value;
        }
    }
    fftw_execute(arrays->q_forward);

    fftw_complex *p_frequency = (fftw_complex *)arrays->p;
    fftw_complex *q_frequency = (fftw_complex *)arrays->q;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < arrays->complex_count; ++i) {
        p_frequency[i][0] = arrays->linear[i] * p_frequency[i][0] +
                            arrays->nonlinear[i] * q_frequency[i][0];
        p_frequency[i][1] = arrays->linear[i] * p_frequency[i][1] +
                            arrays->nonlinear[i] * q_frequency[i][1];
        q_frequency[i][0] = p_frequency[i][0];
        q_frequency[i][1] = p_frequency[i][1];
    }
    fftw_execute(arrays->q_inverse);
}

static void optimize_box(Arrays *arrays, const Model *model,
                         Relaxation *relaxation) {
    double dx0 = relaxation->dx;
    double dy0 = relaxation->dy;
    double d = relaxation->sample_step;
    double dxs[5] = {dx0, dx0 - d, dx0 + d, dx0, dx0};
    double dys[5] = {dy0, dy0, dy0, dy0 - d, dy0 + d};
    double energies[5];
    for (int i = 0; i < 5; ++i) {
        relaxation->dx = dxs[i];
        relaxation->dy = dys[i];
        calculate_properties(arrays, model, relaxation);
        energies[i] = relaxation->free_energy;
    }

    double x_denominator = 2.0 * (energies[1] - 2.0 * energies[0] + energies[2]);
    double y_denominator = 2.0 * (energies[3] - 2.0 * energies[0] + energies[4]);
    relaxation->dx =
        fabs(x_denominator) < 1e-30
            ? dx0
            : (d * (energies[1] - energies[2]) +
               2.0 * dx0 * (-2.0 * energies[0] + energies[1] + energies[2])) /
                  x_denominator;
    relaxation->dy =
        fabs(y_denominator) < 1e-30
            ? dy0
            : (d * (energies[3] - energies[4]) +
               2.0 * dy0 * (-2.0 * energies[0] + energies[3] + energies[4])) /
                  y_denominator;

    double dw = arrays->w * (relaxation->dx - dx0);
    double dh = arrays->h * (relaxation->dy - dy0);
    double change = hypot(dw, dh);
    double limit = PI / sqrt(3.0);
    if (change > limit) {
        double ratio = limit / change;
        relaxation->dx = ratio * relaxation->dx + (1.0 - ratio) * dx0;
        relaxation->dy = ratio * relaxation->dy + (1.0 - ratio) * dy0;
    }
    update_operators(arrays, model, relaxation);

    double dd = hypot(relaxation->dx - dx0, relaxation->dy - dy0);
    relaxation->sample_step *= dd < d ? 0.5 : 2.0;
    if (relaxation->sample_step < 1e-6) {
        relaxation->sample_step *= 2.0;
    }
}

static void write_state(const Arrays *arrays, const Relaxation *relaxation,
                        const Output *output) {
    char filename[320];
    snprintf(filename, sizeof(filename), "%s-t-%d.dat", output->name,
             relaxation->step);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        die("Unable to write state file.");
    }
    int stride = 2 * arrays->wc;
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            fprintf(file, "%.9e\n", arrays->q[(size_t)y * stride + x]);
        }
    }
    fclose(file);
}

static void print_status(const Relaxation *relaxation, const Output *output) {
    printf("%d %lld %.9f %.9f %.9f %.9f\n", relaxation->step,
           (long long)(time(NULL) - relaxation->started), relaxation->dx,
           relaxation->dy, relaxation->free_energy, relaxation->density);
    char filename[320];
    snprintf(filename, sizeof(filename), "%s.out", output->name);
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        die("Unable to write progress output.");
    }
    fprintf(file, "%d %lld %.9f %.9f %.9f %.9f\n", relaxation->step,
            (long long)(time(NULL) - relaxation->started), relaxation->dx,
            relaxation->dy, relaxation->free_energy, relaxation->density);
    fclose(file);
}

static void relax_system(Arrays *arrays, const Model *model,
                         Relaxation *relaxation, const Output *output) {
    for (relaxation->step = 0; relaxation->step <= relaxation->total_steps;
         ++relaxation->step) {
        if (relaxation->optimize_interval > 0 && relaxation->step > 0 &&
            relaxation->step % relaxation->optimize_interval == 0) {
            optimize_box(arrays, model, relaxation);
        }
        if (relaxation->step % output->print_interval == 0) {
            calculate_properties(arrays, model, relaxation);
            print_status(relaxation, output);
        }
        if (relaxation->step % output->write_interval == 0) {
            write_state(arrays, relaxation, output);
        }
        if (relaxation->step < relaxation->total_steps) {
            solver_step(arrays, model, relaxation);
        }
    }
}

static int append_unique(int *ids, int count, int id) {
    for (int i = 0; i < count; ++i) {
        if (ids[i] == id) return count;
    }
    ids[count++] = id;
    return count;
}

static int collect_candidates(const Bin *bins, int nx, int ny, int bx, int by,
                              int *ids) {
    int count = 0;
    for (int oy = -1; oy <= 1; ++oy) {
        int y = wrap_index(by + oy, ny);
        for (int ox = -1; ox <= 1; ++ox) {
            int x = wrap_index(bx + ox, nx);
            const Bin *bin = &bins[y * nx + x];
            for (int i = 0; i < bin->count; ++i) {
                count = append_unique(ids, count, bin->ids[i]);
            }
        }
    }
    return count;
}

static int point_compare(const void *left, const void *right) {
    const Point *a = left;
    const Point *b = right;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    return 0;
}

static void write_xyz_from_field(const Arrays *arrays,
                                 const Relaxation *relaxation,
                                 const char *filename, double pfc_lattice,
                                 double angstrom_lattice) {
    int width = arrays->w;
    int height = arrays->h;
    int stride = 2 * arrays->wc;
    double scale = angstrom_lattice / pfc_lattice;
    double ux = relaxation->dx * scale;
    double uy = relaxation->dy * scale;
    double physical_width = width * ux;
    double physical_height = height * uy;
    double *field = checked_malloc((size_t)width * height * sizeof(*field));
    for (int y = 0; y < height; ++y) {
        memcpy(field + (size_t)y * width,
               arrays->q + (size_t)y * stride, (size_t)width * sizeof(*field));
    }

    PointVector rings = {0};
    const int neighbor_x[8] = {1, 1, 0, -1, -1, -1, 0, 1};
    const int neighbor_y[8] = {0, 1, 1, 1, 0, -1, -1, -1};
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double center = field[(size_t)y * width + x];
            int minimum = 1;
            for (int n = 0; n < 8; ++n) {
                int xx = wrap_index(x + neighbor_x[n], width);
                int yy = wrap_index(y + neighbor_y[n], height);
                if (field[(size_t)yy * width + xx] < center) {
                    minimum = 0;
                    break;
                }
            }
            if (!minimum) continue;

            field[(size_t)y * width + x] -= 1000.0;
            for (int n = 0; n < 8; ++n) {
                int xx = wrap_index(x + neighbor_x[n], width);
                int yy = wrap_index(y + neighbor_y[n], height);
                field[(size_t)yy * width + xx] -= 100.0;
            }

            double p0 = field[(size_t)y * width + wrap_index(x - 1, width)];
            double p1 = field[(size_t)y * width + x];
            double p2 = field[(size_t)y * width + wrap_index(x + 1, width)];
            double denominator = p0 - 2.0 * p1 + p2;
            double refined_x =
                fabs(denominator) < 1e-14
                    ? x * ux
                    : (x + 0.5 * (p0 - p2) / denominator) * ux;

            p0 = field[(size_t)wrap_index(y - 1, height) * width + x];
            p1 = field[(size_t)y * width + x];
            p2 = field[(size_t)wrap_index(y + 1, height) * width + x];
            denominator = p0 - 2.0 * p1 + p2;
            double refined_y =
                fabs(denominator) < 1e-14
                    ? y * uy
                    : (y + 0.5 * (p0 - p2) / denominator) * uy;
            Point ring = {wrap_coordinate(refined_x, physical_width),
                          wrap_coordinate(refined_y, physical_height)};
            point_vector_push(&rings, ring);
        }
    }
    free(field);
    if (rings.count < 3) {
        die("Too few density minima to create carbon coordinates.");
    }

    int nx = (int)ceil(0.25 * physical_width);
    int ny = (int)ceil(0.25 * physical_height);
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    Bin *bins = calloc((size_t)nx * ny, sizeof(*bins));
    if (bins == NULL) die("Out of memory while creating spatial bins.");
    for (size_t i = 0; i < rings.count; ++i) {
        int bx = (int)(rings.items[i].x / physical_width * nx);
        int by = (int)(rings.items[i].y / physical_height * ny);
        if (bx >= nx) bx = nx - 1;
        if (by >= ny) by = ny - 1;
        bin_push(&bins[by * nx + bx], (int)i);
    }

    int max_candidates = 0;
    for (int i = 0; i < nx * ny; ++i) {
        max_candidates += bins[i].count;
    }
    PointVector carbons = {0};
#pragma omp parallel
    {
        PointVector local = {0};
        int *candidates =
            checked_malloc((size_t)max_candidates * sizeof(*candidates));
#pragma omp for schedule(dynamic, 32)
        for (int i = 0; i < (int)rings.count; ++i) {
            Point origin = rings.items[i];
            int bx = (int)(origin.x / physical_width * nx);
            int by = (int)(origin.y / physical_height * ny);
            if (bx >= nx) bx = nx - 1;
            if (by >= ny) by = ny - 1;
            int candidate_count =
                collect_candidates(bins, nx, ny, bx, by, candidates);
            for (int jj = 0; jj < candidate_count; ++jj) {
                int j = candidates[jj];
                if (j <= i) continue;
                Point ij = periodic_difference(rings.items[j], origin,
                                               physical_width, physical_height);
                for (int kk = 0; kk < candidate_count; ++kk) {
                    int k = candidates[kk];
                    if (k <= j) continue;
                    Point ik = periodic_difference(
                        rings.items[k], origin, physical_width, physical_height);
                    double cross = ij.x * ik.y - ij.y * ik.x;
                    if (fabs(cross) < 1e-14) continue;
                    double ij2 = point_norm2(ij);
                    double ik2 = point_norm2(ik);
                    Point circumcenter = {
                        (ij2 * ik.y - ik2 * ij.y) / (2.0 * cross),
                        (ij.x * ik2 - ik.x * ij2) / (2.0 * cross)};
                    double radius2 = point_norm2(circumcenter);
                    int empty = 1;
                    for (int ll = 0; ll < candidate_count; ++ll) {
                        int l = candidates[ll];
                        if (l == i || l == j || l == k) continue;
                        Point il = periodic_difference(
                            rings.items[l], origin, physical_width,
                            physical_height);
                        Point delta = {il.x - circumcenter.x,
                                       il.y - circumcenter.y};
                        if (point_norm2(delta) < radius2 - 1e-12) {
                            empty = 0;
                            break;
                        }
                    }
                    if (empty) {
                        Point carbon = {
                            wrap_coordinate(origin.x + (ij.x + ik.x) / 3.0,
                                            physical_width),
                            wrap_coordinate(origin.y + (ij.y + ik.y) / 3.0,
                                            physical_height)};
                        point_vector_push(&local, carbon);
                    }
                }
            }
        }
#pragma omp critical
        {
            for (size_t i = 0; i < local.count; ++i) {
                point_vector_push(&carbons, local.items[i]);
            }
        }
        free(candidates);
        free(local.items);
    }

    qsort(carbons.items, carbons.count, sizeof(*carbons.items), point_compare);
    FILE *xyz = fopen(filename, "w");
    if (xyz == NULL) {
        die("Unable to write XYZ output.");
    }
    fprintf(xyz, "%zu\n\n", carbons.count);
    for (size_t i = 0; i < carbons.count; ++i) {
        fprintf(xyz, "%zu 1 %.10f %.10f 0\n", i + 1, carbons.items[i].x,
                carbons.items[i].y);
    }
    fclose(xyz);
    printf("Wrote %zu carbon atoms to %s\n", carbons.count, filename);

    for (int i = 0; i < nx * ny; ++i) free(bins[i].ids);
    free(bins);
    free(rings.items);
    free(carbons.items);
}

static void run_input(const char *name, Arrays *arrays, Relaxation *relaxation) {
    char filename[320];
    snprintf(filename, sizeof(filename), "%s.in", name);
    FILE *input = fopen(filename, "r");
    if (input == NULL) {
        fprintf(stderr, "Input file not found: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    Output output = {0};
    Model model = {0};
    memset(arrays, 0, sizeof(*arrays));
    memset(relaxation, 0, sizeof(*relaxation));
    snprintf(output.name, sizeof(output.name), "%s", name);
    relaxation->started = time(NULL);
    relaxation->sample_step = 0.0001;
    snprintf(filename, sizeof(filename), "%s.out", name);
    FILE *progress = fopen(filename, "w");
    if (progress == NULL) die("Unable to create progress output.");
    fclose(progress);

    char label;
    char line[1024];
    while (fscanf(input, " %c", &label) == 1) {
        if (label == '#') {
            (void)fgets(line, sizeof(line), input);
        } else if (label == 'S') {
            seed_rng(input);
        } else if (label == 'O') {
            configure_output(&output, input);
        } else if (label == 'A') {
            configure_arrays(arrays, input);
        } else if (label == 'I') {
            initialize_system(arrays, input);
        } else if (label == 'M') {
            configure_model(&model, input);
        } else if (label == 'R') {
            configure_relaxation(relaxation, input);
            update_operators(arrays, &model, relaxation);
            relax_system(arrays, &model, relaxation, &output);
        } else {
            fprintf(stderr, "Invalid input label: %c\n", label);
            exit(EXIT_FAILURE);
        }
    }
    fclose(input);
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s RUN [RUN ...] [--xyz FILE] [--threads N]\\n"
            "       [--pfc-lattice 7.3] [--angstrom-lattice 2.46]\\n"
            "Example: %s step1 step2 --xyz graphene.xyz --threads 8\\n",
            program, program);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    int run_count = 0;
    int threads = omp_get_max_threads();
    int write_xyz = 1;
    char xyz_filename[320] = "";
    double pfc_lattice = 7.3;
    double angstrom_lattice = 2.46;
    const char **runs = checked_malloc((size_t)argc * sizeof(*runs));
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--xyz") == 0 && i + 1 < argc) {
            snprintf(xyz_filename, sizeof(xyz_filename), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--no-xyz") == 0) {
            write_xyz = 0;
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--pfc-lattice") == 0 && i + 1 < argc) {
            pfc_lattice = atof(argv[++i]);
        } else if (strcmp(argv[i], "--angstrom-lattice") == 0 &&
                   i + 1 < argc) {
            angstrom_lattice = atof(argv[++i]);
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            free(runs);
            return EXIT_FAILURE;
        } else {
            runs[run_count++] = argv[i];
        }
    }
    if (run_count == 0 || threads < 1 || pfc_lattice <= 0.0 ||
        angstrom_lattice <= 0.0) {
        usage(argv[0]);
        free(runs);
        return EXIT_FAILURE;
    }
    if (xyz_filename[0] == '\0') {
        snprintf(xyz_filename, sizeof(xyz_filename), "%s.xyz",
                 runs[run_count - 1]);
    }

    omp_set_num_threads(threads);
    if (!fftw_init_threads()) {
        die("Unable to initialize FFTW threads.");
    }
    fftw_plan_with_nthreads(threads);
    printf("Using %d OpenMP/FFTW threads\n", threads);

    Arrays arrays = {0};
    Relaxation relaxation = {0};
    for (int i = 0; i < run_count; ++i) {
        if (i > 0) clear_arrays(&arrays);
        printf("Running %s\n", runs[i]);
        run_input(runs[i], &arrays, &relaxation);
    }
    if (write_xyz) {
        write_xyz_from_field(&arrays, &relaxation, xyz_filename, pfc_lattice,
                             angstrom_lattice);
    }
    clear_arrays(&arrays);
    fftw_cleanup_threads();
    free(runs);
    return EXIT_SUCCESS;
}
