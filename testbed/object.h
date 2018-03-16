#include "params.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Object {
  public:
    gsl_rng * rng;
    double pos[NDIM];
    void SetRandomPosition() {
      for (int i=0; i<NDIM; ++i) {
        pos[i] = BOX_SIZE*gsl_rng_uniform_pos(rng);
      }
    }
};

