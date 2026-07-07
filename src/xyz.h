#ifndef POLYGR_XYZ_H
#define POLYGR_XYZ_H

typedef enum {
    XYZ_GEOMETRY_PLANE,
    XYZ_GEOMETRY_RIPPLE,
    XYZ_GEOMETRY_BILAYER,
    XYZ_GEOMETRY_TUBE
} XyzGeometry;

typedef struct {
    XyzGeometry geometry;
} XyzOptions;

int xyz_parse_geometry(const char *name, XyzGeometry *geometry);
const char *xyz_geometry_name(XyzGeometry geometry);

/*
 * Convert a row-major PFC density field to graphene carbon coordinates.
 * field_stride is measured in doubles and may be larger than width.
 */
void xyz_write_from_field(const double *field, int width, int height,
                          int field_stride, double dx, double dy,
                          const char *filename, double pfc_lattice,
                          double angstrom_lattice,
                          const XyzOptions *options);

#endif
