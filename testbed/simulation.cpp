#include "simulation.h"

void Simulation::Run() {
  objs.resize(NOBJS);
  InitRNG(); // Init RNGs
  InitPositions(); // Init object positions
}

// Initialize random number generators
void Simulation::InitRNG() {
  gsl_rng_env_setup();
#pragma omp parallel
  {
    int threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    int low = NOBJS*threadnum/numthreads;
    int high = NOBJS*(threadnum+1)/numthreads;
    for (int i=low; i<high; ++i) {
      objs[i].rng = gsl_rng_alloc(gsl_rng_taus);
      gsl_rng_set(objs[i].rng, (i+1)*SEED);
    }
  }
}

void Simulation::InitPositions() {
#pragma omp parallel
  {
    int threadnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    int low = NOBJS*threadnum/numthreads;
    int high = NOBJS*(threadnum+1)/numthreads;
    for (int k=0; k<100; ++k)
    for (int i=low; i<high; ++i) {
      objs[i].SetRandomPosition();;
    }
  }
}

void Simulation::CreatePairs() {

}
