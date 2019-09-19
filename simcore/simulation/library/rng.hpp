#ifndef _SIMCORE_RNG_H_
#define _SIMCORE_RNG_H_

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mutex>
//#include "definitions.hpp"

class RNG {
private:
  const gsl_rng_type *T;
  static long _seed_;
  static std::mutex _rng_mtx_;
  void Clear() { gsl_rng_free(r); }
  void Init() {
    std::lock_guard<std::mutex> lk(_rng_mtx_);
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, _seed_);
    _seed_ = gsl_rng_get(this->r);
  }

public:
  gsl_rng *r;
  RNG() { Init(); }
  ~RNG() { Clear(); }
  static void SetSeed(long seed) { _seed_ = seed; }
  RNG(const RNG &that) : RNG() { 
    gsl_rng_memcpy(this->r, that.r);
  }
  RNG &operator=(RNG const &that) {
    gsl_rng_memcpy(this->r, that.r);
    return *this;
  }
};

void generate_random_unit_vector(int n_dim, double *vec, gsl_rng *r);
// void get_random_coordinate(double * pos, int const n_dim, double radius,
// boundary_type btype, gsl_rng * r);

#endif // _SIMCORE_RNG_H_
