#include "field_image.h"

#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void image_die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(EXIT_FAILURE);
}

static void put_u32_be(FILE *file, uint32_t value) {
    fputc((int)((value >> 24) & 0xff), file);
    fputc((int)((value >> 16) & 0xff), file);
    fputc((int)((value >> 8) & 0xff), file);
    fputc((int)(value & 0xff), file);
}

static uint32_t crc32_update(uint32_t crc, const unsigned char *data,
                             size_t length) {
    crc = ~crc;
    for (size_t i = 0; i < length; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) {
            uint32_t mask = 0u - (crc & 1u);
            crc = (crc >> 1) ^ (0xedb88320u & mask);
        }
    }
    return ~crc;
}

static uint32_t adler32_update(uint32_t adler, const unsigned char *data,
                               size_t length) {
    uint32_t a = adler & 0xffffu;
    uint32_t b = (adler >> 16) & 0xffffu;
    for (size_t i = 0; i < length; ++i) {
        a = (a + data[i]) % 65521u;
        b = (b + a) % 65521u;
    }
    return (b << 16) | a;
}

static void write_chunk(FILE *file, const char type[4],
                        const unsigned char *data, uint32_t length) {
    put_u32_be(file, length);
    fwrite(type, 1, 4, file);
    if (length > 0) {
        fwrite(data, 1, length, file);
    }
    uint32_t crc = crc32_update(0, (const unsigned char *)type, 4);
    crc = crc32_update(crc, data, length);
    put_u32_be(file, crc);
}

static void write_zlib_stored(FILE *file, const unsigned char *data,
                              size_t length) {
    unsigned char *stream =
        malloc(2 + length + 5 * (length / 65535 + 1) + 4);
    if (stream == NULL) image_die("Out of memory while writing PNG.");

    size_t pos = 0;
    stream[pos++] = 0x78;
    stream[pos++] = 0x01;
    size_t remaining = length;
    const unsigned char *cursor = data;
    while (remaining > 0) {
        uint16_t block = remaining > 65535 ? 65535 : (uint16_t)remaining;
        int final = remaining == block;
        stream[pos++] = (unsigned char)(final ? 1 : 0);
        stream[pos++] = (unsigned char)(block & 0xff);
        stream[pos++] = (unsigned char)((block >> 8) & 0xff);
        uint16_t nlen = (uint16_t)~block;
        stream[pos++] = (unsigned char)(nlen & 0xff);
        stream[pos++] = (unsigned char)((nlen >> 8) & 0xff);
        for (uint16_t i = 0; i < block; ++i) {
            stream[pos++] = cursor[i];
        }
        cursor += block;
        remaining -= block;
    }
    uint32_t adler = adler32_update(1, data, length);
    stream[pos++] = (unsigned char)((adler >> 24) & 0xff);
    stream[pos++] = (unsigned char)((adler >> 16) & 0xff);
    stream[pos++] = (unsigned char)((adler >> 8) & 0xff);
    stream[pos++] = (unsigned char)(adler & 0xff);

    write_chunk(file, "IDAT", stream, (uint32_t)pos);
    free(stream);
}

static double sample_normalized(const double *source, int width, int height,
                                int stride, double x, double y,
                                double minimum, double range) {
    if (x < 0.0) x = 0.0;
    if (y < 0.0) y = 0.0;
    if (x > width - 1.0) x = width - 1.0;
    if (y > height - 1.0) y = height - 1.0;
    int x0 = (int)x;
    int y0 = (int)y;
    int x1 = x0 + 1 < width ? x0 + 1 : x0;
    int y1 = y0 + 1 < height ? y0 + 1 : y0;
    double fx = x - x0;
    double fy = y - y0;
    double v00 = source[(size_t)y0 * stride + x0];
    double v10 = source[(size_t)y0 * stride + x1];
    double v01 = source[(size_t)y1 * stride + x0];
    double v11 = source[(size_t)y1 * stride + x1];
    double top = v00 + fx * (v10 - v00);
    double bottom = v01 + fx * (v11 - v01);
    double value = top + fy * (bottom - top);
    double normalized = (value - minimum) / range;
    if (normalized < 0.0) normalized = 0.0;
    if (normalized > 1.0) normalized = 1.0;
    return normalized;
}

void field_write_gray_png(const double *source, int width, int height,
                          int stride, const char *filename, int scale) {
    if (source == NULL || filename == NULL || width <= 0 || height <= 0 ||
        stride < width) {
        image_die("Invalid density field for PNG output.");
    }
    if (scale < 1) scale = 1;
    if (scale > 32) image_die("PNG scale must be between 1 and 32.");

    double minimum = DBL_MAX;
    double maximum = -DBL_MAX;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double value = source[(size_t)y * stride + x];
            if (value < minimum) minimum = value;
            if (value > maximum) maximum = value;
        }
    }
    double range = maximum - minimum;
    if (range == 0.0) range = 1.0;

    int image_width = width * scale;
    int image_height = height * scale;
    size_t row_bytes = (size_t)image_width + 1;
    size_t raw_size = row_bytes * (size_t)image_height;
    unsigned char *raw = malloc(raw_size);
    if (raw == NULL) image_die("Out of memory while rasterizing PNG.");

    for (int y = 0; y < image_height; ++y) {
        unsigned char *row = raw + (size_t)y * row_bytes;
        row[0] = 0;
        double src_y = height - 1.0 - ((double)y + 0.5) / scale;
        for (int x = 0; x < image_width; ++x) {
            double src_x = ((double)x + 0.5) / scale;
            double normalized = sample_normalized(source, width, height,
                                                  stride, src_x, src_y,
                                                  minimum, range);
            row[x + 1] = (unsigned char)(normalized * 255.0 + 0.5);
        }
    }

    FILE *file = fopen(filename, "wb");
    if (file == NULL) image_die("Unable to write PNG output.");
    static const unsigned char signature[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(signature, 1, sizeof(signature), file);

    unsigned char ihdr[13];
    ihdr[0] = (unsigned char)((uint32_t)image_width >> 24);
    ihdr[1] = (unsigned char)((uint32_t)image_width >> 16);
    ihdr[2] = (unsigned char)((uint32_t)image_width >> 8);
    ihdr[3] = (unsigned char)((uint32_t)image_width);
    ihdr[4] = (unsigned char)((uint32_t)image_height >> 24);
    ihdr[5] = (unsigned char)((uint32_t)image_height >> 16);
    ihdr[6] = (unsigned char)((uint32_t)image_height >> 8);
    ihdr[7] = (unsigned char)((uint32_t)image_height);
    ihdr[8] = 8;
    ihdr[9] = 0;
    ihdr[10] = 0;
    ihdr[11] = 0;
    ihdr[12] = 0;
    write_chunk(file, "IHDR", ihdr, sizeof(ihdr));
    write_zlib_stored(file, raw, raw_size);
    write_chunk(file, "IEND", NULL, 0);

    if (fclose(file) != 0) image_die("Unable to finish PNG output.");
    printf("Wrote density PNG to %s\n", filename);
    free(raw);
}

void field_write_rgb_png(const unsigned char *rgb, int width, int height,
                         const char *filename) {
    if (rgb == NULL || filename == NULL || width <= 0 || height <= 0) {
        image_die("Invalid RGB image for PNG output.");
    }

    size_t row_bytes = (size_t)width * 3 + 1;
    size_t raw_size = row_bytes * (size_t)height;
    unsigned char *raw = malloc(raw_size);
    if (raw == NULL) image_die("Out of memory while rasterizing RGB PNG.");

    for (int y = 0; y < height; ++y) {
        unsigned char *row = raw + (size_t)y * row_bytes;
        row[0] = 0;
        memcpy(row + 1, rgb + (size_t)y * width * 3, (size_t)width * 3);
    }

    FILE *file = fopen(filename, "wb");
    if (file == NULL) image_die("Unable to write RGB PNG output.");
    static const unsigned char signature[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(signature, 1, sizeof(signature), file);

    unsigned char ihdr[13];
    ihdr[0] = (unsigned char)((uint32_t)width >> 24);
    ihdr[1] = (unsigned char)((uint32_t)width >> 16);
    ihdr[2] = (unsigned char)((uint32_t)width >> 8);
    ihdr[3] = (unsigned char)((uint32_t)width);
    ihdr[4] = (unsigned char)((uint32_t)height >> 24);
    ihdr[5] = (unsigned char)((uint32_t)height >> 16);
    ihdr[6] = (unsigned char)((uint32_t)height >> 8);
    ihdr[7] = (unsigned char)((uint32_t)height);
    ihdr[8] = 8;
    ihdr[9] = 2;
    ihdr[10] = 0;
    ihdr[11] = 0;
    ihdr[12] = 0;
    write_chunk(file, "IHDR", ihdr, sizeof(ihdr));
    write_zlib_stored(file, raw, raw_size);
    write_chunk(file, "IEND", NULL, 0);

    if (fclose(file) != 0) image_die("Unable to finish RGB PNG output.");
    printf("Wrote figure PNG to %s\n", filename);
    free(raw);
}
