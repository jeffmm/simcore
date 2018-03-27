#include <iostream>
#include "src/rng.h"
#include <vector>

#define BOX_SIZE 100
#define OBJ_NUMBER 10
#define SEED 419694816
#define NDIM 2


class Obj {
  private:
    double pos_[NDIM];
    static long seed_;
    RNG rng_;
  public:
    Obj() {
      rng_.Init(seed_)
      seed_ = 1000*gsl_rng_get(rng_.r)-gsl_rng_get(rng_.r);//gsl_rng_get(rng_.r);
      for (int i=0; i<NDIM; ++i) {
        pos_[i] = BOX_SIZE*gsl_rng_uniform_pos(rng_.r);
      }
    }

};
long Obj::seed_ = SEED;

int main() {
  std::vector<Obj> objs;
  objs.resize(OBJ_NUMBER);
  //CreatePairs();
  return 0;
}
