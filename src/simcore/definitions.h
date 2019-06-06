#ifndef _SIMCORE_DEFINITIONS_H_
#define _SIMCORE_DEFINITIONS_H_

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
  public:                                      \
    Enum() = default;
#include "enum.h"

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

BETTER_ENUM(species_id, unsigned char, br_bead, filament, passive_filament, centrosome, bead_spring, spherocylinder, spindle, motor)
BETTER_ENUM(draw_type, unsigned char, fixed, orientation, bw, none);
BETTER_ENUM(potential_type, unsigned char, none, wca, soft);
BETTER_ENUM(boundary_type, unsigned char, none=0, box=1, sphere=2, budding=3);
BETTER_ENUM(poly_state, unsigned char, grow, shrink, pause);

struct space_struct {
  int n_dim;
  int n_periodic;
  bool bud;
  double radius;
  double pressure_tensor[9]; // pressure tensor
  double pressure; // isometric pressure
  double volume;
  double bud_radius;
  double bud_height;
  double bud_neck_radius;
  double bud_neck_height;
  double *unit_cell;
  double *unit_cell_inv; // inverse unit cell
  double *a; // direct lattice vector
  double *b; // reciprocal lattice vector
  double *a_perp; // perpendicular distance between opposite unit cell faces
  double *mu; // scaling matrix for constant pressure
  int n_bound; // number of bound motors
  double concentration; // C of motors
  boundary_type type;
};

struct graph_struct {
  double r[3];
  double u[3];
  double length;
  double diameter;
  double color;
  draw_type draw;
};

#endif
