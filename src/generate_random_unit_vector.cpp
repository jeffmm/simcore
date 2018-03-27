/* This function takes either a 2d or 3d vector and places it randomly on the surface
   of an n_dim sphere */

#include <gsl/gsl_rng.h>
#include <math.h>

void generate_random_unit_vector(int n_dim, double *vect, gsl_rng *r) {
    double x, y, z, w, t;

    w = 1.0;
    if (n_dim == 3) {
        z = 2.0 * gsl_rng_uniform_pos(r) - 1.0;
        w = sqrt(1 - z * z);
        vect[2] = z;
    }

    t = 2.0 * M_PI * gsl_rng_uniform_pos(r);
    x = w * cos(t);
    y = w * sin(t);
    vect[0] = x;
    vect[1] = y;
}

