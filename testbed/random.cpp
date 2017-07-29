#include "rng.h"

void generate_random_unit_vector(double * vec, int const n_dim, gsl_rng * r) {
    double w = 1.0;
    if (n_dim == 3) {
      double z = 2.0 * gsl_rng_uniform_pos(r) - 1.0;
      w = sqrt(1-z*z);
      vec[2] = z;
    }
    double t = 2.0 * M_PI * gsl_rng_uniform_pos(r);
    double x = w * cos(t);
    double y = w * sin(t);
    vec[0] = x;
    vec[1] = y;
}

void get_random_coordinate(double * pos, int const n_dim, double radius, boundary_type btype, gsl_rng * r) {
  if (btype == BOX) {
    for (int i=0; i<n_dim; ++i) {
      pos[i] = radius * (2.0*gsl_rng_uniform_pos(r)-1);
    }
  }
  else if (btype == SPHERE) {
    generate_random_unit_vector(pos, n_dim, r);
    for (int i=0; i<n_dim; ++i) {
      pos[i] *= radius;
    }
  }
}


