#include "xyz.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static void xyz_die(const char *message) {
    fprintf(stderr, "%s\n", message);
    abort();
}

static void *xyz_malloc(size_t bytes) {
    void *ptr = malloc(bytes);
    if (ptr == NULL) xyz_die("Out of memory during XYZ coordinate extraction.");
    return ptr;
}

static void point_vector_push(PointVector *vector, Point point) {
    if (vector->count == vector->capacity) {
        size_t next = vector->capacity == 0 ? 1024 : vector->capacity * 2;
        Point *items = realloc(vector->items, next * sizeof(*items));
        if (items == NULL) {
            xyz_die("Out of memory while storing carbon coordinates.");
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
        if (ids == NULL) xyz_die("Out of memory while building spatial bins.");
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

static Point periodic_difference(Point b, Point a, double width,
                                 double height) {
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

static double point_norm2(Point point) {
    return point.x * point.x + point.y * point.y;
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

void xyz_write_from_field(const double *source, int width, int height,
                          int field_stride, double dx, double dy,
                          const char *filename, double pfc_lattice,
                          double angstrom_lattice,
                          const XyzOptions *options) {
    (void)options;
    if (source == NULL || filename == NULL || width < 2 || height < 2 ||
        field_stride < width || dx <= 0.0 || dy <= 0.0 ||
        pfc_lattice <= 0.0 || angstrom_lattice <= 0.0) {
        xyz_die("Invalid density field or lattice parameters for XYZ output.");
    }

    double scale = angstrom_lattice / pfc_lattice;
    double ux = dx * scale;
    double uy = dy * scale;
    double physical_width = width * ux;
    double physical_height = height * uy;
    double *field = xyz_malloc((size_t)width * height * sizeof(*field));
    for (int y = 0; y < height; ++y) {
        memcpy(field + (size_t)y * width,
               source + (size_t)y * field_stride,
               (size_t)width * sizeof(*field));
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
        xyz_die("Too few density minima to create carbon coordinates.");
    }

    int nx = (int)ceil(0.25 * physical_width);
    int ny = (int)ceil(0.25 * physical_height);
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    Bin *bins = calloc((size_t)nx * ny, sizeof(*bins));
    if (bins == NULL) xyz_die("Out of memory while creating spatial bins.");
    for (size_t i = 0; i < rings.count; ++i) {
        int bx = (int)(rings.items[i].x / physical_width * nx);
        int by = (int)(rings.items[i].y / physical_height * ny);
        if (bx >= nx) bx = nx - 1;
        if (by >= ny) by = ny - 1;
        bin_push(&bins[by * nx + bx], (int)i);
    }

    int max_candidates = (int)rings.count;
    PointVector carbons = {0};
#pragma omp parallel
    {
        PointVector local = {0};
        int *candidates =
            xyz_malloc((size_t)max_candidates * sizeof(*candidates));
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
    if (xyz == NULL) xyz_die("Unable to write XYZ output.");
    fprintf(xyz, "%zu\n", carbons.count);
    fprintf(xyz,
            "Lattice=\"%.10f 0 0 0 %.10f 0 0 0 20\" "
            "Origin=\"0 0 -10\" "
            "Properties=\"id:I:1:species:S:1:pos:R:3:"
            "radius:R:1:color:R:3\" "
            "pbc=\"T T F\"\n",
            physical_width, physical_height);
    for (size_t i = 0; i < carbons.count; ++i) {
        fprintf(xyz, "%zu C %.10f %.10f 0 0.70 0.35 0.35 0.35\n", i + 1,
                carbons.items[i].x, carbons.items[i].y);
    }
    fclose(xyz);
    printf("Wrote %zu carbon atoms to %s\n", carbons.count, filename);

    for (int i = 0; i < nx * ny; ++i) free(bins[i].ids);
    free(bins);
    free(rings.items);
    free(carbons.items);
}
