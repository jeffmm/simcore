#ifndef _SIMCORE_RNG_H_
#define _SIMCORE_RNG_H_

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "auxiliary.hpp"

class RNG {
 private:
  gsl_rng *rng_;

 public:
  RNG(unsigned long seed);
  ~RNG();
  RNG(const RNG &that) : RNG(that.GetSeed()) {}
  RNG &operator=(RNG const &that) { return *this; }
  const double RandomUniform();
  const int RandomPoisson(const double mean);
  const long RandomInt(const long n);
  const double RandomNormal(const double sigma);
  void RandomUnitVector(const int n_dim, double *vec);
  void RandomCoordinate(const space_struct *const s, double *vec,
                        const double buffer = 0);
  unsigned long GetSeed() const;
  void *GetState();
  size_t GetSize();
  template <typename T>
  void Shuffle(T *array, size_t size);
};

template <typename T>
void RNG::Shuffle(T *array, size_t size) {
  gsl_ran_shuffle(rng_, array, size, sizeof(T));
}

#endif  // _SIMCORE_RNG_H_
