#include "field_image.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double x;
    double y;
    int gb;
} Atom;

typedef struct {
    Atom *items;
    size_t count;
    size_t capacity;
    double width;
    double height;
} XyzData;

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t bytes) {
    void *ptr = malloc(bytes);
    if (ptr == NULL) die("Out of memory.");
    return ptr;
}

static void push_atom(XyzData *data, Atom atom) {
    if (data->count == data->capacity) {
        size_t next = data->capacity == 0 ? 1024 : data->capacity * 2;
        Atom *items = realloc(data->items, next * sizeof(*items));
        if (items == NULL) die("Out of memory while reading XYZ atoms.");
        data->items = items;
        data->capacity = next;
    }
    data->items[data->count++] = atom;
}

static void parse_lattice(const char *line, double *width, double *height) {
    const char *p = strstr(line, "Lattice=\"");
    if (p == NULL) die("Extended XYZ header is missing Lattice.");
    p += strlen("Lattice=\"");
    double values[9];
    if (sscanf(p, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &values[0], &values[1], &values[2], &values[3], &values[4],
               &values[5], &values[6], &values[7], &values[8]) != 9) {
        die("Unable to parse Extended XYZ lattice.");
    }
    *width = values[0];
    *height = values[4];
    if (*width <= 0.0 || *height <= 0.0) {
        die("Invalid Extended XYZ lattice dimensions.");
    }
}

static void read_xyz(const char *path, XyzData *data) {
    FILE *file = fopen(path, "r");
    if (file == NULL) die("Unable to open XYZ input.");

    char line[4096];
    if (fgets(line, sizeof(line), file) == NULL) die("XYZ is missing atom count.");
    long atom_count = strtol(line, NULL, 10);
    if (atom_count <= 0) die("XYZ atom count is not positive.");
    if (fgets(line, sizeof(line), file) == NULL) die("XYZ is missing header.");
    parse_lattice(line, &data->width, &data->height);

    for (long i = 0; i < atom_count; ++i) {
        if (fgets(line, sizeof(line), file) == NULL) die("XYZ ended early.");
        int id, coordination, ring0, ring1, ring2, gb;
        char species[16];
        double x, y, z, radius, red, green, blue;
        int fields = sscanf(line,
                            "%d %15s %lf %lf %lf %d %d %d %d %d %lf %lf %lf %lf",
                            &id, species, &x, &y, &z, &coordination, &ring0,
                            &ring1, &ring2, &gb, &radius, &red, &green, &blue);
        (void)id;
        (void)species;
        (void)z;
        (void)coordination;
        (void)ring0;
        (void)ring1;
        (void)ring2;
        (void)radius;
        (void)red;
        (void)green;
        (void)blue;
        if (fields != 14) die("XYZ atom row does not match gb-enabled format.");
        Atom atom = {x, y, gb};
        push_atom(data, atom);
    }
    fclose(file);
}

static void pixel(unsigned char *rgb, int width, int height, int x, int y,
                  int r, int g, int b) {
    if (x < 0 || y < 0 || x >= width || y >= height) return;
    unsigned char *p = rgb + ((size_t)y * width + x) * 3;
    p[0] = (unsigned char)r;
    p[1] = (unsigned char)g;
    p[2] = (unsigned char)b;
}

static void disc(unsigned char *rgb, int width, int height, int cx, int cy,
                 int radius, int r, int g, int b) {
    for (int y = cy - radius; y <= cy + radius; ++y) {
        for (int x = cx - radius; x <= cx + radius; ++x) {
            int dx = x - cx;
            int dy = y - cy;
            if (dx * dx + dy * dy <= radius * radius) {
                pixel(rgb, width, height, x, y, r, g, b);
            }
        }
    }
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

static double periodic_delta(double value, double box) {
    if (value > 0.5 * box) value -= box;
    if (value < -0.5 * box) value += box;
    return value;
}

static void map_xy(const XyzData *data, int image_w, int image_h, int margin,
                   double x, double y, int *px, int *py) {
    double sx = (image_w - 2.0 * margin) / data->width;
    double sy = (image_h - 2.0 * margin) / data->height;
    *px = margin + (int)floor(x * sx + 0.5);
    *py = image_h - margin - (int)floor(y * sy + 0.5);
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

static void write_gb_map(const XyzData *data, const char *output, int image_w,
                         int image_h, int grid_x, int grid_y,
                         double connect_cutoff) {
    int margin = 28;
    unsigned char *rgb = xmalloc((size_t)image_w * image_h * 3);
    memset(rgb, 255, (size_t)image_w * image_h * 3);

    if (grid_x > 0 && grid_y > 0) {
        draw_grid(rgb, image_w, image_h, margin, grid_x, grid_y);
    }

    double cutoff = connect_cutoff;
    double cutoff2 = cutoff * cutoff;
    for (size_t i = 0; i < data->count; ++i) {
        if (!data->items[i].gb) continue;
        int x0, y0;
        map_xy(data, image_w, image_h, margin, data->items[i].x,
               data->items[i].y, &x0, &y0);
        for (size_t j = i + 1; j < data->count; ++j) {
            if (!data->items[j].gb) continue;
            double dx = periodic_delta(data->items[j].x - data->items[i].x,
                                       data->width);
            double dy = periodic_delta(data->items[j].y - data->items[i].y,
                                       data->height);
            if (dx * dx + dy * dy > cutoff2) continue;
            double xj = data->items[i].x + dx;
            double yj = data->items[i].y + dy;
            if (xj < 0.0 || xj >= data->width || yj < 0.0 ||
                yj >= data->height) {
                continue;
            }
            int x1, y1;
            map_xy(data, image_w, image_h, margin, xj, yj, &x1, &y1);
            line(rgb, image_w, image_h, x0, y0, x1, y1, 238, 72, 61);
        }
    }

    int gb_count = 0;
    for (size_t i = 0; i < data->count; ++i) {
        if (!data->items[i].gb) continue;
        ++gb_count;
        int x, y;
        map_xy(data, image_w, image_h, margin, data->items[i].x,
               data->items[i].y, &x, &y);
        disc(rgb, image_w, image_h, x, y, 1, 238, 72, 61);
    }

    field_write_rgb_png(rgb, image_w, image_h, output);
    printf("Wrote grain-boundary map with %d gb atoms to %s\n", gb_count,
           output);
    free(rgb);
}

static void usage(const char *program) {
    fprintf(stderr,
            "Usage: %s input.xyz --output gb-map.png [--width 512] "
            "[--height 512] [--grid NX NY] [--connect CUTOFF_ANGSTROM]\n",
            program);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    const char *input = NULL;
    const char *output = "gb-map.png";
    int width = 512;
    int height = 512;
    int grid_x = 4;
    int grid_y = 4;
    double connect_cutoff = 1.95;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output = argv[++i];
        } else if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            width = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--height") == 0 && i + 1 < argc) {
            height = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--grid") == 0 && i + 2 < argc) {
            grid_x = atoi(argv[++i]);
            grid_y = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--connect") == 0 && i + 1 < argc) {
            connect_cutoff = atof(argv[++i]);
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
    if (input == NULL || width < 64 || height < 64 || grid_x < 0 ||
        grid_y < 0 || connect_cutoff <= 0.0) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    XyzData data = {0};
    read_xyz(input, &data);
    write_gb_map(&data, output, width, height, grid_x, grid_y, connect_cutoff);
    free(data.items);
    return EXIT_SUCCESS;
}
