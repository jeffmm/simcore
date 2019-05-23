#ifndef _SIMCORE_RNG_H_
#define _SIMCORE_RNG_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
//#include "definitions.h"

class RNG { 
  private:
    const gsl_rng_type *T;
    void Clear() {
      gsl_rng_free(r);
    }
  public:
    gsl_rng *r;
    RNG() {}
    RNG(long seed) {
      Init(seed);
    }
    ~RNG() {
      Clear();
    }
    void Init(long seed) {
      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);
      gsl_rng_set(r, seed);
    }
    RNG(const RNG& that) {
      this->Init(gsl_rng_get(that.r));
    }
    RNG& operator=(RNG const&that) {
      this->Init(gsl_rng_get(that.r));
      return *this;
    }
};

void generate_random_unit_vector(int n_dim, double * vec, gsl_rng * r);
//void get_random_coordinate(double * pos, int const n_dim, double radius, boundary_type btype, gsl_rng * r);

#endif // _SIMCORE_RNG_H_
