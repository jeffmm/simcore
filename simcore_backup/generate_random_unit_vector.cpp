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

/* This function takes two unit vectors and rotates vect1 such that its relative
 * z axis is aligned with unit vector vect2. This allows the placement of an
 * orientation vector on the surface of a sphere such that it is orientated away
 * from the sphere (ie we wont have any randomly placed mts on a centrosome that
 * are pointing inward). The formulas for rotation was derived using rodrigues */

void rotate_orientation_vector(int n_dim, double *vect1, double *vect2) {
  double vx, vy, vz, theta, phi;
  vx = vect1[0];
  vy = vect1[1];
  if (n_dim == 3) vz = vect1[2];

  // get theta, phi for vect2
  phi = atan2(vect2[1],vect2[0]);
  if (n_dim == 3) theta = acos(vect2[2]);

  // rodrigues formula
  if (n_dim == 3) {
    vect1[0] = vx * cos(theta) * cos(phi) - vy * sin(phi) + vz * cos(phi) * sin(theta);
    vect1[1] = vx * cos(theta) * sin(phi) + vy * cos(phi) + vz * sin(theta) * sin(phi);
    vect1[2] = vz * cos(theta) - vx * sin(theta);
  }
  else if (n_dim==2) {
    vect1[0] = vx * cos(phi) - vy * sin(phi);
    vect1[1] = vx * sin(phi) + vy * cos(phi);
  }
}

