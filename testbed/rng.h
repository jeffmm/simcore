#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

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

