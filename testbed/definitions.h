#ifndef _SIMCORE_DEFINITIONS_H_
#define _SIMCORE_DEFINITIONS_H_

#include <string>

enum boundary_type {
  BOX,
  SPHERE 
};

struct graph_struct {
  double r[3];
  double u[3]; 
  double color[4];
  double length;
  double diameter;
  int draw_type;  // 0: use color array, 1: color by orientation
};

// Structure that holds all relevant spatial data
struct space_struct {
  int n_dim;
  int n_periodic;
  double radius; // system radius
  double *unit_cell;
  double *unit_cell_inv; // inverse unit cell
  double *a; // direct lattice vector
  double *b; // reciprocal lattice vector
  double *a_perp; // perpendicular distance between opposite unit cell faces
  // Info for drawing a budding yeast boundary
  bool bud;
  double bud_radius;
  double bud_height;
  std::string type; // boundary type (sphere, box, or budding)
};

#endif
