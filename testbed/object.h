#include "params.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class Object {
  public:
    gsl_rng * rng;
    double pos[NDIM];
    std::vector<int> nlist;
    void SetRandomPosition() {
      for (int i=0; i<NDIM; ++i) {
        pos[i] = BOX_SIZE*(gsl_rng_uniform_pos(rng)-0.5);
      }
    }
    Object() {
      //for (int i=0;i<NDIM;++i) {
        //pos[i] = 0.0;
      //}
      //nlist.reserve(50);
    }
};

