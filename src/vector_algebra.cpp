#include <math.h>
#include "auxiliary.h"

double dot_product(int n_dim, double *a, double *b)
{
    int i;
    double mag = 0.0;
    for (i = 0; i < n_dim; ++i)
        mag += a[i] * b[i];

    return mag;
}

void cross_product(double *a, double *b, double *c, int n_dim) {
  if (n_dim == 2) {
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
  else if (n_dim == 3) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
}

void normalize_vector(double *a, int n_dim) {
  double mag = 0.0;
  for (int i=0; i<n_dim; ++i) 
    mag += a[i]*a[i];
  mag = sqrt(mag);
  for (int i=0; i<n_dim; ++i) 
    a[i]/=mag;
}
 
// Rotates vector a counterclockwise around a UNIT VECTOR b,
// by an amount theta
void rotate_3d_vector(double theta, double *a, double *b) {
  double t[3];
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  for (int i=0; i<3; ++i)
    t[i] = a[i];
  a[0] = t[0]*(b[0]*b[0]+(b[1]*b[1]+b[2]*b[2])*cos_theta)
        +t[2]*(b[0]*b[2]*(1-cos_theta)+b[1]*sin_theta)
        +t[1]*(b[0]*b[1]*(1-cos_theta)-b[2]*sin_theta);
  a[1] = t[1]*(b[1]*b[1]+(b[0]*b[0]+b[2]*b[2])*cos_theta)
        +t[2]*(b[1]*b[2]*(1-cos_theta)-b[0]*sin_theta)
        +t[0]*(b[0]*b[1]*(1-cos_theta)+b[2]*sin_theta);
  a[2] = t[2]*(b[2]*b[2]+(b[1]*b[1]+b[0]*b[0])*cos_theta)
        +t[1]*(b[1]*b[2]*(1-cos_theta)+b[0]*sin_theta)
        +t[0]*(b[0]*b[2]*(1-cos_theta)-b[1]*sin_theta); 
}

//Inverts symmetric "linear" 2D matrix
//Input:  2x2 Matrix to be inverted (a)
//Output: Inverted matrix ( b)
void invert_sym_2d_matrix(double *a, double *b) {
    double inv_det;
    inv_det = a[0]*a[3] - a[1]*a[1];
    inv_det = 1.0 / inv_det;
    b[0] = a[3] * inv_det;
    b[3] = a[0] * inv_det;
    b[1] = b[2] = -a[1] * inv_det;

    return;
}

//Inverts symmetric "linear" 3D matrix
//Input:  3x3 Matrix to be inverted (a)
//Output: Inverted matrix ( b)
void invert_sym_3d_matrix(double *a, double *b) {
    double inv_det;
    inv_det = a[0] * (a[4]*a[8] - a[5]*a[7])
        - a[1] * (a[3]*a[8] - a[5]*a[6])
        + a[2] * (a[3]*a[7] - a[4]*a[6]);
    inv_det = 1.0/inv_det;
    b[0] = inv_det * (a[4]*a[8] - a[5]*a[7]);
    b[4] = inv_det * (a[0]*a[8] - a[2]*a[6]);
    b[8] = inv_det * (a[0]*a[4] - a[1]*a[3]);
    b[1] = b[3] = -inv_det * (a[3]*a[8] - a[5]*a[6]);
    b[2] = b[6] = inv_det * (a[3]*a[7] - a[4]*a[6]);
    b[5] = b[7] = -inv_det * (a[0]*a[7] - a[1]*a[6]);

    return;
}

