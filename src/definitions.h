#ifndef _SIMCORE_DEFINITIONS_H_
#define _SIMCORE_DEFINITIONS_H_

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
  public:                                      \
    Enum() = default;
#include "enum.h"

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

BETTER_ENUM(species_id, unsigned char, br_bead, filament, centrosome)
BETTER_ENUM(draw_type, unsigned char, fixed, orientation, bw, none);
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

// Harmonic anchor structure
class Anchor {
  public:
    bool alignment_potential_;
    double position_[3],
           orientation_[3],
           force_[3],
           torque_[3],
           k_spring_,
           k_align_,
           spring_length_;
    Anchor() {
      std::fill(position_,position_+3,0.0);
      std::fill(orientation_,orientation_+3,0.0);
      std::fill(force_,force_+3,0.0);
      std::fill(torque_,torque_+3,0.0);
      alignment_potential_ = false;
      k_spring_ = 0;
      k_align_ = 0;
      spring_length_ = 0;
    }
    void ZeroForce() {
      std::fill(force_,force_+3,0.0);
      std::fill(torque_,torque_+3,0.0);
    }
};

typedef std::vector<Anchor>::iterator anchor_iterator;

#endif
