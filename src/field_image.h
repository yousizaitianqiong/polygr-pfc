#ifndef FIELD_IMAGE_H
#define FIELD_IMAGE_H

void field_write_gray_png(const double *source, int width, int height,
                          int stride, const char *filename, int scale);
void field_write_rgb_png(const unsigned char *rgb, int width, int height,
                         const char *filename);

#endif
