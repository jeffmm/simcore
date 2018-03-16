#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include "omp.h"
#include "object.h"


int ix(int x, int y);
void print_array(double * a);
void init_rng(std::vector<Object> * objs);
void populate(std::vector<Object> * objs);

int main() {
  std::vector<Object> objs;
  objs.resize(NOBJS);
  init_rng(&objs); // Init RNGs
  populate(&objs); // Init object positions

  // **** DO STUFF HERE **** //

  // **** DO STUFF HERE **** //

  return 0;
}

// Initialize random number generators
void init_rng(std::vector<Object> * objs) {
  gsl_rng_env_setup();
#pragma omp parallel
  {
    int threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    int low = NOBJS*threadnum/numthreads;
    int high = NOBJS*(threadnum+1)/numthreads;
    for (int i=low; i<high; ++i) {
      (*objs)[i].rng = gsl_rng_alloc(gsl_rng_taus);
      gsl_rng_set((*objs)[i].rng, (i+1)*SEED);
    }
  }
}

void populate(std::vector<Object> * objs) {
#pragma omp parallel
  {
    int threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    int low = NOBJS*threadnum/numthreads;
    int high = NOBJS*(threadnum+1)/numthreads;
    for (int k=0; k<100; ++k)
    for (int i=low; i<high; ++i) {
      (*objs)[i].SetRandomPosition();;
    }
  }
}

void print_array(double * a) {
  printf("{ ");
  for (int i=0; i<NOBJS*NDIM; ++i) {
    printf("%1.3f ",a[i]);
  }
  printf("}\n");
}

int ix(int x, int y) {
  return x*NDIM + y;
}

