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

#include "xyz.h"

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
            char *comment = fgets(line, sizeof(line), input);
            (void)comment;
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
            "Usage: %s RUN [RUN ...] [--xyz FILE] [--threads N]\n"
            "       [--pfc-lattice 7.3] [--angstrom-lattice 2.46]\n"
            "Example: %s step1 step2 --xyz graphene.xyz --threads 8\n",
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
        xyz_write_from_field(arrays.q, arrays.w, arrays.h, 2 * arrays.wc,
                             relaxation.dx, relaxation.dy, xyz_filename,
                             pfc_lattice, angstrom_lattice);
    }
    clear_arrays(&arrays);
    fftw_cleanup_threads();
    free(runs);
    return EXIT_SUCCESS;
}
