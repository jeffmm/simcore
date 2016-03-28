#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <gsl/gsl_rng.h>
#include <vector>
#include <iostream>
#include "allocate.h"
#include "macros.h"
#include "timetester.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "math.h"
#include "parameters.h"
#include <ios>
#include <fstream>
#include <sstream>

// No u

struct graph_struct {
  int n_spheros;
  int i_sphero;
  double **h;
  double **r;
  double **u;
  double *l;
  double *diam;
  double m_rad;
  double d_rad;
  double m_d_dist;
};

struct rng_properties { 
  gsl_rng *r;
  const gsl_rng_type *T;
  void init(long seed) {
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
  }
  void clear() {
    gsl_rng_free(r);
  }
};

struct interaction {
  int tid; // type of interacting object (to determine potential)
  double *dr; // distance to object
  double rad; // interaction radius of object
};

double cpu();
void generate_random_unit_vector(int n_dim, double *vect, gsl_rng *r);
void rotate_orientation_vector(int n_dim, double *vect1, double *vect2);
double dot_product(int n_dim, double *a, double *b);
double *separation_vector(int n_dim, int n_periodic, double *r1, double *s1, double *r2, double *s2, double **unit_cell);
void cross_product(double *a, double *b, double *c, int n_dim);
void normalize_vector(double *a, int n_dim);
void error_exit(const char *error_msg, ...);
void warning(const char *warning_msg);
void tridiagonal_solver(double *a, double *b, double *c, double *d, int n);
void rotate_3d_vector(double theta, double *a, double *b);
void invert_sym_2d_matrix(double *a, double *b); 
void invert_sym_3d_matrix(double *a, double *b); 

#endif // _AUXILIARY_H_

