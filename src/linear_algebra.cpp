#include <simcore/auxiliary.hpp>

double dot_product(int n_dim, double const *const a, double const *const b) {
  double mag = 0.0;
  for (int i = 0; i < n_dim; ++i) {
    mag += a[i] * b[i];
  }
  return mag;
}

// WARNING: Must always use n_dim=3 if you are using torques,
//          since they point in the z direction when using 2D
void cross_product(double const *const a, double const *const b, double *c,
                   int n_dim) {
  if (n_dim == 2) {
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = a[0] * b[1] - a[1] * b[0];
  } else if (n_dim == 3) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
}

void normalize_vector(double *a, int n_dim) {
  double mag = 0.0;
  for (int i = 0; i < n_dim; ++i) {
    mag += a[i] * a[i];
  }
  mag = sqrt(mag);
  for (int i = 0; i < n_dim; ++i) {
    a[i] /= mag;
  }
}

// Inverts symmetric "linear" 2D matrix
// Input:  2x2 Matrix to be inverted (a)
// Output: Inverted matrix ( b)
void invert_sym_2d_matrix(double *a, double *b) {
  double inv_det;
  inv_det = a[0] * a[3] - a[1] * a[1];
  inv_det = 1.0 / inv_det;
  b[0] = a[3] * inv_det;
  b[3] = a[0] * inv_det;
  b[2] = b[1] = -a[1] * inv_det;
}

// Inverts symmetric "linear" 3D matrix
// Input:  3x3 Matrix to be inverted (a)
// Output: Inverted matrix ( b)
void invert_sym_3d_matrix(double *a, double *b) {
  double inv_det;
  inv_det = a[0] * (a[4] * a[8] - a[5] * a[7]) -
            a[1] * (a[3] * a[8] - a[5] * a[6]) +
            a[2] * (a[3] * a[7] - a[4] * a[6]);
  inv_det = 1.0 / inv_det;
  b[0] = inv_det * (a[4] * a[8] - a[5] * a[7]);
  b[4] = inv_det * (a[0] * a[8] - a[2] * a[6]);
  b[8] = inv_det * (a[0] * a[4] - a[1] * a[3]);
  b[1] = b[3] = -inv_det * (a[3] * a[8] - a[5] * a[6]);
  b[2] = b[6] = inv_det * (a[3] * a[7] - a[4] * a[6]);
  b[5] = b[7] = -inv_det * (a[0] * a[7] - a[1] * a[6]);
}

void separation_vector(int n_dim, int n_periodic, double const *const r1,
                       double const *const s1, double const *const r2,
                       double const *const s2, double *unit_cell, double *dr) {
  // First handle periodic subspace
  double ds[3];
  for (int i = 0; i < n_periodic; ++i) {
    ds[i] = s2[i] - s1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (int i = 0; i < n_periodic; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j) {
      dr[i] += unit_cell[n_dim * i + j] * ds[j];
    }
  }
  // Then handle free subspace
  for (int i = n_periodic; i < n_dim; ++i) {
    dr[i] = r2[i] - r1[i];
  }
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
void periodic_boundary_conditions(int n_dim, int n_periodic, double *h,
                                  double *h_inv, double *r, double *s) {
  /* Compute scaled coordinate and apply periodic boundary conditions. */
  for (int i = 0; i < n_periodic; ++i) {
    s[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j) {
      s[i] += h_inv[n_dim * i + j] * r[j];
    }
    s[i] -= NINT(s[i]);
  }

  /* Recompute real coordinates accounting for periodic boundary conditions. */
  // for (int i = 0; i < n_periodic; ++i) {
  // r[i] = 0.0;
  // for (int j = 0; j < n_periodic; ++j)
  // r[i] += h[n_dim*i+j] * s[j];
  //}
}

/* Computes determinant of square matrix mat of size n using recursion (and
   ugly dynamic memory allocation ) */
double determinant(int n, double **mat) {
  double det = 0;
  double **submat = new double *[n];
  for (int i = 0; i < n; ++i) {
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
    for (int c = 0; c < n; ++c) {
      int subi = 0;
      for (int i = 1; i < n; ++i) {
        int subj = 0;
        for (int j = 0; j < n; ++j) {
          if (j == c) {
            continue;
          }
          submat[subi][subj] = mat[i][j];
          subj++;
        }
        subi++;
      }
      det = det + (pow(-1, c) * mat[0][c] * determinant(n - 1, submat));
    }
  }
  for (int i = 0; i < n; ++i) {
    delete[] submat[i];
  }
  delete[] submat;
  return det;
}

void tridiagonal_solver(std::vector<double> *a, std::vector<double> *b,
                        std::vector<double> *c, std::vector<double> *d, int n) {
  //  Solves Ax=d for x as seen below:
  //
  //  | b[0] c[0]   0  | | x[0] |   | d[0] |
  //  | a[0] b[1] c[1] | | x[1] | = | d[1] |
  //  |   0  a[1] b[2] | | x[2] |   | d[2] |

  //  Works only for tridiagonal matrices.
  //  Uses Thomas' algorithm which is unstable
  //  unless matrix is guaranteed to be symmetric
  //  positive definite.

  //  Input: pointer to array of lower off-diagonal elements (*a)
  //         pointer to array of diagonal elements as seen above (*b)
  //         pointer to array of upper off-diagonal elements  (*c)
  //         pointer to array of d vector (*d)
  //         number of unknowns (n)
  //
  //  Output: Array of solutions x (*d)

  n--;
  (*c)[0] /= (*b)[0];
  (*d)[0] /= (*b)[0];
  for (int i = 1; i < n; ++i) {
    (*c)[i] /= (*b)[i] - (*a)[i - 1] * (*c)[i - 1];
    (*d)[i] = ((*d)[i] - (*a)[i - 1] * (*d)[i - 1]) /
              ((*b)[i] - (*a)[i - 1] * (*c)[i - 1]);
  }
  (*d)[n] = ((*d)[n] - (*a)[n - 1] * (*d)[n - 1]) /
            ((*b)[n] - (*a)[n - 1] * (*c)[n - 1]);
  for (int i = n; i-- > 0;) (*d)[i] -= (*c)[i] * (*d)[i + 1];
  return;
}

/* This function rotates vector v about vector k by an angle theta.
 * Derived using rodrigues' rotation formula */
void rotate_vector(double *v, double *k, double theta, int n_dim) {
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  if (n_dim == 2) {
    v[0] = v[0] * cos_theta - v[1] * sin_theta;
    v[1] = v[0] * sin_theta + v[1] * cos_theta;
  } else {
    double k_dot_v = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    double t[3];
    t[0] = v[0] * cos_theta + (k[1] * v[2] - v[1] * k[2]) * sin_theta +
           k[0] * k_dot_v * (1 - cos_theta);
    t[1] = v[1] * cos_theta + (k[2] * v[0] - k[0] * v[2]) * sin_theta +
           k[1] * k_dot_v * (1 - cos_theta);
    t[2] = v[2] * cos_theta + (k[0] * v[1] - k[1] * v[0]) * sin_theta +
           k[2] * k_dot_v * (1 - cos_theta);

    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[2];
  }
}

/* This function takes two unit vectors and rotates vect1 such that its relative
 * z axis is aligned with unit vector vect2. This allows the placement of an
 * orientation vector on the surface of a sphere such that it is orientated away
 * from the sphere (ie we wont have any randomly placed mts on a centrosome that
 * are pointing inward). The formulas for rotation was derived using rodrigues
 */

void rotate_vector_relative(int n_dim, double *vect1, double *vect2) {
  double vx, vy, vz, theta, phi;
  vx = vect1[0];
  vy = vect1[1];
  if (n_dim == 3) vz = vect1[2];

  // get theta, phi for vect2
  phi = atan2(vect2[1], vect2[0]);
  if (n_dim == 3) theta = acos(vect2[2]);

  // rodrigues formula
  if (n_dim == 3) {
    vect1[0] =
        vx * cos(theta) * cos(phi) - vy * sin(phi) + vz * cos(phi) * sin(theta);
    vect1[1] =
        vx * cos(theta) * sin(phi) + vy * cos(phi) + vz * sin(theta) * sin(phi);
    vect1[2] = vz * cos(theta) - vx * sin(theta);
  } else if (n_dim == 2) {
    vect1[0] = vx * cos(phi) - vy * sin(phi);
    vect1[1] = vx * sin(phi) + vy * cos(phi);
  }
}
