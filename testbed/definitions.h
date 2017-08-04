#ifndef _SIMCORE_DEFINITIONS_H_
#define _SIMCORE_DEFINITIONS_H_

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string>

#define SQR(x)             ((x) * (x))
#define CUBE(x)            ((x) * (x) * (x))
#define NINT(x)            ((x) < 0.0 ? (int) ((x) - 0.5) : (int) ((x) + 0.5))
#define ABS(x)             ((x) < 0 ? -(x) : (x))
#define MAX(x,y)           ((x) > (y) ? (x) : (y))
#define MIN(x,y)           ((x) < (y) ? (x) : (y))
#define SIGN(a,b)          ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SIGNOF(x)          ((x) >= 0.0 ? 1 : -1)


enum boundary_type {
  BOX,
  SPHERE 
};

struct graph_struct {
  double r[3];
  double u[3]; 
  double color;
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

double dot_product(int n_dim, double const * const a, double const * const b);
void cross_product(double const * const a, double const * const b, double *c, int n_dim);
void normalize_vector(double *a, int n_dim);

#endif
