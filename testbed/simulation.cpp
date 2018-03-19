#include "simulation.h"

void Simulation::Run() {
  objs.resize(NOBJS);
  InitRNG(); // Init RNGs
  InitPositions(); // Init object positions
  CreatePairs();
  //PrintPairs();
}

// Initialize random number generators
void Simulation::InitRNG() {
  gsl_rng_env_setup();
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      objs[i].rng = gsl_rng_alloc(gsl_rng_taus);
      gsl_rng_set(objs[i].rng, (i+1)*SEED);
    }
  }
}

void Simulation::InitPositions() {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    //for (int k=0; k<100; ++k) // For time testing
    for (int i=low; i<high; ++i) {
      objs[i].SetRandomPosition();;
    }
  }
}

//void Simulation::CreatePairs() {

//}
