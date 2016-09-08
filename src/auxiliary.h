#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <array>
#include <map>
#include <tuple>
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

#include "helpers.h"

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

extern bool debug_trace;

struct graph_struct {
  double r[3];
  double u[3];
  double length;
  double diameter;
};

enum class SID : unsigned char {
  none,
  md_bead,
  md_bead_opt,
  br_bead,
  br_dimer,
  argon,
  neon,
  br_rod,
  filament,
  br_simple_rod,
  md_kmc_bead,
  br_walker,
  xlink_head,
  xlink
};

//XXX: CJE I HATE THIS WITH THE FUCKING DYING PASSION
//OF A THOUSAND SUNS!!!!
inline SID StringToSID(std::string &s) {
    if (s == "none")
        return SID::none;
    else if (s == "md_bead")
        return SID::md_bead;
    else if (s == "md_bead_opt")
        return SID::md_bead_opt;
    else if (s == "br_bead")
        return SID::br_bead;
    else if (s == "br_dimer")
        return SID::br_dimer;
    else if (s == "argon")
        return SID::argon;
    else if (s == "neon")
        return SID::neon;
    else if (s == "br_rod")
        return SID::br_rod;
    else if (s == "filament")
        return SID::filament;
    else if (s == "br_simple_rod")
        return SID::br_simple_rod;
    else if (s == "md_kmc_bead")
        return SID::md_kmc_bead;
    else if (s == "br_walker")
        return SID::br_walker;
    else if (s == "xlink_head")
        return SID::xlink_head;
    else if (s == "xlink")
        return SID::xlink;
    else
        std::cout << "Go fuck yourself!\n";
    return SID::none;
}

//XXX: And again but backwards and even more annoying
inline std::string SIDToString(SID &s) {
    if (s == SID::none)
        return "none";
    else if (s == SID::md_bead)
        return "md_bead";
    else if (s == SID::md_bead_opt)
        return "md_bead_opt";
    else if (s == SID::br_bead)
        return "br_bead";
    else if (s == SID::br_dimer)
        return "br_dimer";
    else if (s == SID::argon)
        return "argon";
    else if (s == SID::neon)
        return "neon";
    else if (s == SID::br_rod)
        return "br_rod";
    else if (s == SID::filament)
        return "filament";
    else if (s == SID::br_simple_rod)
        return "br_simple_rod";
    else if (s == SID::md_kmc_bead)
        return "md_kmc_bead";
    else if (s == SID::br_walker)
        return "br_walker";
    else if (s == SID::xlink_head)
        return "xlink_head";
    else if (s == SID::xlink)
        return "xlink";
    else
        std::cout << "Go fuck yourself!\n";
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

struct space_struct {
  int n_dim;
  int n_periodic;
  bool bud;
  double radius;
  double bud_radius;
  double bud_height;
  double **unit_cell;
  double **unit_cell_inv;
  double **a;
  double **b;
  double *a_perp;
  std::string type;
};


double cpu();
void generate_random_unit_vector(int n_dim, double *vect, gsl_rng *r);
void dummy_function();
void rotate_orientation_vector(int n_dim, double *vect1, double *vect2);
double dot_product(int n_dim, double const * const a, double const * const b);
double determinant(int n, double **mat);
//double *separation_vector(int n_dim, int n_periodic, double *r1, double *s1, double *r2, double *s2, double **unit_cell);
void separation_vector(int n_dim, int n_periodic, double const * const r1, double const * const s1, double const * const r2, double const * const s2, double **unit_cell, double *dr);
void cross_product(double *a, double *b, double *c, int n_dim);
void normalize_vector(double *a, int n_dim);
void error_exit(const char *error_msg, ...);
void warning(const char *warning_msg);
void tridiagonal_solver(std::vector<double> *a, std::vector<double> *b, std::vector<double> *c, std::vector<double> *d, int n);
void rotate_3d_vector(double theta, double *a, double *b);
void invert_sym_2d_matrix(double **a, double **b); 
void invert_sym_3d_matrix(double **a, double **b); 
void periodic_boundary_conditions(int n_periodic, double **h, double **h_inv,
                                         double *r, double *s);

class PotentialBase;
struct interaction {
  PotentialBase *potential; // interaction potential (returns force)
  double dr[3]; // separation vector from this obj to other obj
  double contact[3]; // separation vector from position_ to point of contact
  double buffer; // sum of 'skin depth' of both objects (sum of radii for spheres)
  double dr_mag; //magnitude of dr separation vector
};

// XXX: CJE working title for now
struct interactionmindist {
    double dr_mag;
    double dr_mag2;
    double buffer_mag;
    double buffer_mag2;
    double dr[3];
    double contact1[3];
    double contact2[3];
};

namespace cytohelpers {
    // From bithacks online
    __attribute__((always_inline))
    static inline unsigned int nextpow2(unsigned int v) {
        v--;

        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;

        return v;
    }

    // Cell vector to linear id
    __attribute__((always_inline))
    static inline int cell_vec_to_linear(int cx, int cy, int cz, int nc[3]) {
        return cx + cy*nc[0] + cz*nc[0]*nc[1];
    }

    // Cell linear id to vector id
    __attribute__((always_inline))
    static inline void cell_linear_to_vec(int cidx, int nc[3], int* cx) {
        cx[0] = cidx % nc[0];
        cx[1] = (cidx / nc[0]) % nc[1];
        cx[2] = cidx / (nc[0] * nc[1]);
    }

}

#endif // _AUXILIARY_H_

