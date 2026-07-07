/*
 * Single-GPU CUDA/cuFFT implementation of the shared-memory PFC solver path.
 *
 * This backend intentionally keeps the time-step loop on the GPU. Host copies
 * are used for initialization, progress/output files, and box optimization.
 */

#include <cuda_runtime.h>
#include <cufft.h>

#include <math.h>
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
    int stride;
    size_t real_count;
    size_t complex_count;
    double *host_q;
    double *host_p;
    double *linear;
    double *nonlinear;
    double *p;
    double *q;
    cufftHandle p_forward;
    cufftHandle q_forward;
    cufftHandle q_inverse;
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

static void check_cuda(cudaError_t status, const char *where) {
    if (status != cudaSuccess) {
        fprintf(stderr, "%s: %s\n", where, cudaGetErrorString(status));
        exit(EXIT_FAILURE);
    }
}

static void check_cufft(cufftResult status, const char *where) {
    if (status != CUFFT_SUCCESS) {
        fprintf(stderr, "%s: cuFFT status %d\n", where, (int)status);
        exit(EXIT_FAILURE);
    }
}

static void *checked_malloc(size_t bytes) {
    void *ptr = malloc(bytes);
    if (ptr == NULL) die("Out of memory.");
    return ptr;
}

__global__ static void update_operators_kernel(
    double *linear_out, double *nonlinear_out, int w, int h, int wc,
    double alpha, double beta, double conserved, double dx, double dy,
    double dt) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t count = (size_t)h * wc;
    if (i >= count) return;
    int y = (int)(i / wc);
    int x = (int)(i - (size_t)y * wc);
    double dkx = 2.0 * PI / (dx * w);
    double dky = 2.0 * PI / (dy * h);
    double kx = x * dkx;
    double ky = (y < h / 2 ? y : y - h) * dky;
    double k2 = kx * kx + ky * ky;
    double one_minus_k2 = 1.0 - k2;
    double linear = alpha + beta * one_minus_k2 * one_minus_k2;
    double mobility = (1.0 - conserved) + conserved * k2;
    double exponential = exp(-mobility * linear * dt);
    double inv_size = 1.0 / ((double)w * h);
    linear_out[i] = exponential;
    nonlinear_out[i] =
        (linear == 0.0 ? -mobility * dt : (exponential - 1.0) / linear) *
        inv_size;
}

__global__ static void scale_complex_kernel(cufftDoubleComplex *data,
                                            size_t count, double scale) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    data[i].x *= scale;
    data[i].y *= scale;
}

__global__ static void nonlinear_kernel(double *q, int w, int h, int stride,
                                        double gamma, double delta) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t count = (size_t)w * h;
    if (i >= count) return;
    int y = (int)(i / w);
    int x = (int)(i - (size_t)y * w);
    size_t index = (size_t)y * stride + x;
    double value = q[index];
    q[index] = gamma * value * value + delta * value * value * value;
}

__global__ static void frequency_update_kernel(
    cufftDoubleComplex *p, cufftDoubleComplex *q, const double *linear,
    const double *nonlinear, size_t count) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    double real = linear[i] * p[i].x + nonlinear[i] * q[i].x;
    double imag = linear[i] * p[i].y + nonlinear[i] * q[i].y;
    p[i].x = real;
    p[i].y = imag;
    q[i].x = real;
    q[i].y = imag;
}

__global__ static void spectral_factor_kernel(cufftDoubleComplex *q, int w,
                                              int h, int wc, double dx,
                                              double dy) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t count = (size_t)h * wc;
    if (i >= count) return;
    int y = (int)(i / wc);
    int x = (int)(i - (size_t)y * wc);
    double dkx = 2.0 * PI / (dx * w);
    double dky = 2.0 * PI / (dy * h);
    double kx = x * dkx;
    double ky = (y < h / 2 ? y : y - h) * dky;
    double factor = 1.0 - kx * kx - ky * ky;
    factor *= factor;
    q[i].x *= factor;
    q[i].y *= factor;
}

static int blocks_for(size_t count) {
    return (int)((count + 255) / 256);
}

static void configure_arrays(Arrays *arrays, FILE *input) {
    if (fscanf(input, " %d %d", &arrays->w, &arrays->h) != 2 ||
        arrays->w < 2 || arrays->h < 2) {
        die("Invalid array dimensions.");
    }
    arrays->wc = arrays->w / 2 + 1;
    arrays->stride = 2 * arrays->wc;
    arrays->real_count = (size_t)arrays->h * arrays->stride;
    arrays->complex_count = (size_t)arrays->h * arrays->wc;
    arrays->host_q = (double *)checked_malloc(arrays->real_count * sizeof(double));
    arrays->host_p = (double *)checked_malloc(arrays->real_count * sizeof(double));
    memset(arrays->host_q, 0, arrays->real_count * sizeof(double));
    memset(arrays->host_p, 0, arrays->real_count * sizeof(double));
    check_cuda(cudaMalloc(&arrays->linear, arrays->complex_count * sizeof(double)), "cudaMalloc linear");
    check_cuda(cudaMalloc(&arrays->nonlinear, arrays->complex_count * sizeof(double)), "cudaMalloc nonlinear");
    check_cuda(cudaMalloc(&arrays->p, arrays->real_count * sizeof(double)), "cudaMalloc p");
    check_cuda(cudaMalloc(&arrays->q, arrays->real_count * sizeof(double)), "cudaMalloc q");
    check_cufft(cufftPlan2d(&arrays->p_forward, arrays->h, arrays->w, CUFFT_D2Z), "cufftPlan2d p_forward");
    check_cufft(cufftPlan2d(&arrays->q_forward, arrays->h, arrays->w, CUFFT_D2Z), "cufftPlan2d q_forward");
    check_cufft(cufftPlan2d(&arrays->q_inverse, arrays->h, arrays->w, CUFFT_Z2D), "cufftPlan2d q_inverse");
}

static void clear_arrays(Arrays *arrays) {
    if (arrays->p_forward) cufftDestroy(arrays->p_forward);
    if (arrays->q_forward) cufftDestroy(arrays->q_forward);
    if (arrays->q_inverse) cufftDestroy(arrays->q_inverse);
    cudaFree(arrays->p);
    cudaFree(arrays->q);
    cudaFree(arrays->linear);
    cudaFree(arrays->nonlinear);
    free(arrays->host_q);
    free(arrays->host_p);
    memset(arrays, 0, sizeof(*arrays));
}

static void seed_rng(FILE *input) {
    unsigned int seed;
    if (fscanf(input, " %u", &seed) != 1) die("Invalid random seed.");
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
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            arrays->host_q[(size_t)y * arrays->stride + x] =
                density + amplitude * ((double)rand() / RAND_MAX - 0.5);
        }
    }
}

static void polycrystalline_state(Arrays *arrays, double dx, double dy,
                                  double lattice, double density,
                                  double amplitude, int grain_count,
                                  double radius) {
    if (grain_count <= 0 || radius <= 0.0) die("Invalid polycrystalline initialization.");
    double *grains = (double *)checked_malloc((size_t)grain_count * 3 * sizeof(double));
    for (int grain = 0; grain < grain_count; ++grain) {
        grains[3 * grain] = (double)rand() / RAND_MAX * arrays->w;
        grains[3 * grain + 1] = (double)rand() / RAND_MAX * arrays->h;
        grains[3 * grain + 2] = 2.0 * PI * rand() / RAND_MAX;
    }
    double divisor = 0.5 / (radius * radius * lattice * lattice);
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
            arrays->host_q[(size_t)y * arrays->stride + x] =
                density + amplitude * exp(-closest_r2 * divisor) *
                              one_mode(closest_x, closest_y, lattice,
                                       grains[3 * closest + 2]);
        }
    }
    free(grains);
}

static void read_state(Arrays *arrays, FILE *file, double density,
                       double amplitude) {
    double mean = 0.0;
    double minimum = HUGE_VAL;
    double maximum = -HUGE_VAL;
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            double value;
            if (fscanf(file, "%lf", &value) != 1) die("Invalid state data.");
            arrays->host_q[(size_t)y * arrays->stride + x] = value;
            mean += value;
            if (value < minimum) minimum = value;
            if (value > maximum) maximum = value;
        }
    }
    mean /= (double)arrays->w * arrays->h;
    double range = maximum - minimum;
    if (range == 0.0) die("Cannot normalize a constant input state.");
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            size_t index = (size_t)y * arrays->stride + x;
            arrays->host_q[index] =
                (arrays->host_q[index] - mean) / range * amplitude + density;
        }
    }
}

static void initialize_system(Arrays *arrays, FILE *input) {
    int type;
    if (fscanf(input, " %d", &type) != 1) die("Invalid initialization type.");
    if (type == 0) {
        double density, amplitude;
        if (fscanf(input, " %lf %lf", &density, &amplitude) != 2) die("Invalid random initialization.");
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
        if (fscanf(input, " %255s %lf %lf", filename, &density, &amplitude) != 3) {
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
    check_cuda(cudaMemcpy(arrays->q, arrays->host_q,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyHostToDevice), "copy initial q");
}

static void update_operators(Arrays *arrays, const Model *model,
                             const Relaxation *relaxation) {
    update_operators_kernel<<<blocks_for(arrays->complex_count), 256>>>(
        arrays->linear, arrays->nonlinear, arrays->w, arrays->h, arrays->wc,
        model->alpha, model->beta, model->conserved, relaxation->dx,
        relaxation->dy, relaxation->dt);
    check_cuda(cudaGetLastError(), "update_operators_kernel");
}

static void scale_complex(double *data, size_t count, double scale) {
    scale_complex_kernel<<<blocks_for(count), 256>>>((cufftDoubleComplex *)data,
                                                     count, scale);
    check_cuda(cudaGetLastError(), "scale_complex_kernel");
}

static void calculate_properties(Arrays *arrays, const Model *model,
                                 Relaxation *relaxation) {
    check_cuda(cudaMemcpy(arrays->p, arrays->q,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyDeviceToDevice), "copy q to p");
    check_cufft(cufftExecD2Z(arrays->q_forward, arrays->q,
                             (cufftDoubleComplex *)arrays->q), "q_forward");
    scale_complex(arrays->q, arrays->complex_count,
                  1.0 / ((double)arrays->w * arrays->h));
    spectral_factor_kernel<<<blocks_for(arrays->complex_count), 256>>>(
        (cufftDoubleComplex *)arrays->q, arrays->w, arrays->h, arrays->wc,
        relaxation->dx, relaxation->dy);
    check_cuda(cudaGetLastError(), "spectral_factor_kernel");
    check_cufft(cufftExecZ2D(arrays->q_inverse,
                             (cufftDoubleComplex *)arrays->q, arrays->q),
                "q_inverse properties");
    check_cuda(cudaMemcpy(arrays->host_p, arrays->p,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyDeviceToHost), "copy p host");
    check_cuda(cudaMemcpy(arrays->host_q, arrays->q,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyDeviceToHost), "copy q host");

    double free_energy = 0.0;
    double density = 0.0;
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            size_t index = (size_t)y * arrays->stride + x;
            double p = arrays->host_p[index];
            double p2 = p * p;
            free_energy +=
                0.5 * (model->alpha * p2 +
                       model->beta * p * arrays->host_q[index]) +
                model->gamma * p * p2 / 3.0 +
                0.25 * model->delta * p2 * p2;
            density += p;
            arrays->host_q[index] = p;
        }
    }
    double inv_size = 1.0 / ((double)arrays->w * arrays->h);
    relaxation->free_energy = free_energy * inv_size;
    relaxation->density = density * inv_size;
    check_cuda(cudaMemcpy(arrays->q, arrays->host_q,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyHostToDevice), "restore q");
    check_cufft(cufftExecD2Z(arrays->p_forward, arrays->p,
                             (cufftDoubleComplex *)arrays->p), "p_forward");
    scale_complex(arrays->p, arrays->complex_count, inv_size);
}

static void solver_step(Arrays *arrays, const Model *model,
                        const Relaxation *relaxation) {
    if (relaxation->step % 100 == 0) {
        check_cuda(cudaMemcpy(arrays->p, arrays->q,
                              arrays->real_count * sizeof(double),
                              cudaMemcpyDeviceToDevice), "solver copy q to p");
        check_cufft(cufftExecD2Z(arrays->p_forward, arrays->p,
                                 (cufftDoubleComplex *)arrays->p), "solver p_forward");
        scale_complex(arrays->p, arrays->complex_count,
                      1.0 / ((double)arrays->w * arrays->h));
    }
    nonlinear_kernel<<<blocks_for((size_t)arrays->w * arrays->h), 256>>>(
        arrays->q, arrays->w, arrays->h, arrays->stride, model->gamma,
        model->delta);
    check_cuda(cudaGetLastError(), "nonlinear_kernel");
    check_cufft(cufftExecD2Z(arrays->q_forward, arrays->q,
                             (cufftDoubleComplex *)arrays->q), "solver q_forward");
    frequency_update_kernel<<<blocks_for(arrays->complex_count), 256>>>(
        (cufftDoubleComplex *)arrays->p, (cufftDoubleComplex *)arrays->q,
        arrays->linear, arrays->nonlinear, arrays->complex_count);
    check_cuda(cudaGetLastError(), "frequency_update_kernel");
    check_cufft(cufftExecZ2D(arrays->q_inverse,
                             (cufftDoubleComplex *)arrays->q, arrays->q),
                "solver q_inverse");
}

static void write_state(Arrays *arrays, const Relaxation *relaxation,
                        const Output *output) {
    check_cuda(cudaMemcpy(arrays->host_q, arrays->q,
                          arrays->real_count * sizeof(double),
                          cudaMemcpyDeviceToHost), "copy state host");
    char filename[320];
    snprintf(filename, sizeof(filename), "%s-t-%d.dat", output->name,
             relaxation->step);
    FILE *file = fopen(filename, "w");
    if (file == NULL) die("Unable to write state file.");
    for (int y = 0; y < arrays->h; ++y) {
        for (int x = 0; x < arrays->w; ++x) {
            fprintf(file, "%.9e\n", arrays->host_q[(size_t)y * arrays->stride + x]);
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
    if (file == NULL) die("Unable to write progress output.");
    fprintf(file, "%d %lld %.9f %.9f %.9f %.9f\n", relaxation->step,
            (long long)(time(NULL) - relaxation->started), relaxation->dx,
            relaxation->dy, relaxation->free_energy, relaxation->density);
    fclose(file);
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
    if (relaxation->sample_step < 1e-6) relaxation->sample_step *= 2.0;
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
            "Usage: %s RUN [RUN ...] [--no-xyz] [--no-png]\n"
            "CUDA backend currently writes .out and .dat state files only.\n",
            program);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    int device = 0;
    int write_xyz = 1;
    int write_png = 1;
    int run_count = 0;
    const char **runs = (const char **)checked_malloc((size_t)argc * sizeof(char *));
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--device") == 0 && i + 1 < argc) {
            device = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--no-xyz") == 0) {
            write_xyz = 0;
        } else if (strcmp(argv[i], "--no-png") == 0) {
            write_png = 0;
        } else if (strcmp(argv[i], "--xyz") == 0 || strcmp(argv[i], "--png") == 0 ||
                   strcmp(argv[i], "--threads") == 0 ||
                   strcmp(argv[i], "--png-scale") == 0 ||
                   strcmp(argv[i], "--pfc-lattice") == 0 ||
                   strcmp(argv[i], "--angstrom-lattice") == 0) {
            ++i;
        } else if (argv[i][0] == '-') {
            usage(argv[0]);
            free(runs);
            return EXIT_FAILURE;
        } else {
            runs[run_count++] = argv[i];
        }
    }
    if (run_count == 0) {
        usage(argv[0]);
        free(runs);
        return EXIT_FAILURE;
    }
    if (write_xyz || write_png) {
        fprintf(stderr, "polygr_cuda currently supports solver validation only; pass --no-xyz --no-png.\n");
        free(runs);
        return EXIT_FAILURE;
    }
    check_cuda(cudaSetDevice(device), "cudaSetDevice");
    cudaDeviceProp prop;
    check_cuda(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
    printf("Using CUDA device %d: %s\n", device, prop.name);

    Arrays arrays = {0};
    Relaxation relaxation = {0};
    for (int i = 0; i < run_count; ++i) {
        if (i > 0) clear_arrays(&arrays);
        printf("Running %s\n", runs[i]);
        run_input(runs[i], &arrays, &relaxation);
    }
    clear_arrays(&arrays);
    free(runs);
    return EXIT_SUCCESS;
}
