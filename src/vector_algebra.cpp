#include <math.h>
#include "auxiliary.h"

double dot_product(int n_dim, double const * const a, double const * const b)
{
    int i;
    double mag = 0.0;
    for (i = 0; i < n_dim; ++i)
        mag += a[i] * b[i];

    return mag;
}

// WARNING: Must always use n_dim=3 if you are using torques,
//          since they point in the z direction when using 2D
void cross_product(double const * const a, double const * const b, double *c, int n_dim) {
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
    b[2] = b[1] = -a[1] * inv_det;
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
}

void separation_vector(int n_dim, int n_periodic, double const * const r1, double const * const s1, double const * const r2, double const * const s2, double *unit_cell, double *dr) {

  // First handle periodic subspace
  double ds[3];
  for (int i = 0; i < n_periodic; ++i) {  
    ds[i] = s2[i] - s1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (int i = 0; i < n_periodic; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j)
      dr[i] += unit_cell[n_dim*i+j] * ds[j];
  }
  // Then handle free subspace
  for (int i = n_periodic; i < n_dim; ++i) 
      dr[i] = r2[i] - r1[i];
  return;
}

/* 
Inputs:
   n_dim, number of dimensions
   n_periodic, number of periodic dimensions
   h, unit cell matrix
   h_inv, inverse unit cell matrix
   r, position
   s, scaled position
Outputs:
   new position
   new scaled position
   */
void periodic_boundary_conditions(int n_dim, int n_periodic, double *h, double *h_inv,
                                         double *r, double *s) {
  /* Compute scaled coordinate and apply periodic boundary conditions. */
  for (int i = 0; i < n_periodic; ++i) {
    s[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j) 
      s[i] += h_inv[n_dim*i+j] * r[j];
    s[i] -= NINT(s[i]);
  }

  /* Recompute real coordinates accounting for periodic boundary conditions. */
  for (int i = 0; i < n_periodic; ++i) {
    r[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j) 
      r[i] += h[n_dim*i+j] * s[j];
  }
}

/* Computes determinant of square matrix mat of size n using recursion (and
   ugly dynamic memory allocation ) */
double determinant(int n, double **mat) {
  double det=0;
  double **submat = new double*[n];
  for (int i=0; i<n; ++i) {
    submat[i] = new double[n];
  }
  if (n == 1) {
    for (int i = 0; i < n; ++i) {
      delete[] submat[i];
    }
    delete[] submat;
    return mat[0][0];
  }
  if (n == 2) {
    for (int i = 0; i < n; ++i) {
      delete[] submat[i];
    }
    delete[] submat;
    return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
  } else {
    for(int c = 0; c < n; ++c) {
      int subi = 0;
      for(int i = 1; i < n; ++i) {
        int subj = 0;
        for(int j = 0; j < n; ++j) {
          if (j == c) {
            continue;
          }
          submat[subi][subj] = mat[i][j];
          subj++;
        }
        subi++;
      }
      det = det + (pow(-1 ,c) * mat[0][c] * determinant(n - 1 ,submat));
    }
  }
  for (int i=0; i<n; ++i)
    delete[] submat[i];
  delete[] submat;
  return det;
}

void dummy_function() {return;}
