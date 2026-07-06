#ifndef POLYGR_XYZ_H
#define POLYGR_XYZ_H

/*
 * Convert a row-major PFC density field to graphene carbon coordinates.
 * field_stride is measured in doubles and may be larger than width.
 */
void xyz_write_from_field(const double *field, int width, int height,
                          int field_stride, double dx, double dy,
                          const char *filename, double pfc_lattice,
                          double angstrom_lattice);

#endif
