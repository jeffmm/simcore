#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <array>
#include <map>
#include <tuple>
#include <iostream>
#include <iomanip>
#include "allocate.h"
#include "macros.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "math.h"
#include "parameters.h"
#include <ios>
#include <fstream>
#include <sstream>

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

extern bool debug_trace;
extern bool early_exit;

struct graph_struct {
  double r[3];
  double u[3];
  double color[4];
  double length;
  double diameter;
  int draw_type;
};

//NOTE: Species ids cannot be more than 20 characters for reasons -AL
enum class SID : unsigned char {
  none,
  md_bead,
  br_bead,
  hard_rod,
  filament,
};

inline SID StringToSID(std::string &s) {
  if (s == "none")
    return SID::none;
  else if (s == "md_bead")
    return SID::md_bead;
  else if (s == "br_bead")
    return SID::br_bead;
  else if (s == "hard_rod")
    return SID::hard_rod;
  else if (s == "filament")
    return SID::filament;
  else {
    printf("ERROR: SID string %s not recognized\n",s.c_str());
    exit(1);
  }
  return SID::none;
}

inline std::string SIDToString(SID const s) {
  if (s == SID::none)
    return "none";
  else if (s == SID::md_bead)
    return "md_bead";
  else if (s == SID::br_bead)
    return "br_bead";
  else if (s == SID::hard_rod)
    return "hard_rod";
  else if (s == SID::filament)
    return "filament";
  else {
    printf("SID value not recognized\n");
    exit(1);
  }
  return "none";
}

enum class FTYPE : unsigned char {
  none,
  allpairs,
  microcells,
  cells,
  neighborallpairs,
  neighborcells
};

typedef std::pair<SID, SID> sid_pair;

class rng_properties { 
  private:
    const gsl_rng_type *T;
    void clear() {
      gsl_rng_free(r);
    }
  public:
    gsl_rng *r;
    rng_properties() {}
    rng_properties(long seed) {
      init(seed);
    }
    ~rng_properties() {
      clear();
    }
    void init(long seed) {
      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);
      gsl_rng_set(r, seed);
    }
    rng_properties(const rng_properties& that) {
      this->init(gsl_rng_get(that.r));
    }
    rng_properties& operator=(rng_properties const&that) {
      this->init(gsl_rng_get(that.r));
      return *this;
    }
};

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
  std::string type;
};


double cpu();
void grabber(int width, int height, std::string fname, int framenum);
void generate_random_unit_vector(int n_dim, double *vect, gsl_rng *r);
void rotate_orientation_vector(int n_dim, double *vect1, double *vect2);
double dot_product(int n_dim, double const * const a, double const * const b);
double determinant(int n, double **mat);
void separation_vector(int n_dim, int n_periodic, double const * const r1, double const * const s1, double const * const r2, double const * const s2, double **unit_cell, double *dr);
void cross_product(double const * const a, double const * const b, double *c, int n_dim);
void normalize_vector(double *a, int n_dim);
void error_exit(const char *error_msg, ...);
void warning(const char *warning_msg);
void tridiagonal_solver(std::vector<double> *a, std::vector<double> *b, std::vector<double> *c, std::vector<double> *d, int n);
void rotate_3d_vector(double theta, double *a, double *b);
void invert_sym_2d_matrix(double *a, double *b); 
void invert_sym_3d_matrix(double *a, double *b); 
void periodic_boundary_conditions(int n_dim, int n_periodic, double *h, double *h_inv,
                                         double *r, double *s);

#endif // _AUXILIARY_H_

