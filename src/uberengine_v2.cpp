#include <chrono>

#include "object.h"
#include "particle_engine.h"
#include "uberengine_v2.h"

// Pass in the main system properties information
void UberEngineV2::Init(system_parameters *pParams,
                        space_struct *pSpace,
                        std::vector<SpeciesBase*> *pSpecies,
                        al_set *pAnchors,
                        long seed) {
  if(debug_trace) {
    std::cout << "UberEngineV2::Init\n";
  }
  params_ = pParams;
  space_ = pSpace;
  species_ = pSpecies;
  anchors_ = pAnchors;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  rng_.init(seed);
  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  {
    if (0 == omp_get_thread_num()) {
      nthreads_ = omp_get_num_threads();
    }
  }
  #else
  nthreads_ = 1;
  #endif

  // Create the particle engine scheme, which will take care of the potentials as well
  pengine_.Init(params_, space_, species_, anchors_, &interactions_, gsl_rng_get(rng_.r));
  pengine_.CreateTracking();
  //pengine_.UpdateInteractions();

  fengine_.Init(space_, &pengine_, &interactions_);
  fengine_.InitMP();

  kengine_.Init(params_, space_, &pengine_, &interactions_, gsl_rng_get(rng_.r));
  kengine_.InitMP();

  pengine_.UpdateInteractions();

  fengine_.Interact();

  pengine_.Print();
  kengine_.Print();
  pengine_.Dump();
  fengine_.Dump();
  kengine_.Dump();
}

void UberEngineV2::InteractMP() {
  pengine_.UpdateInteractions();
  fengine_.Interact();
}

void UberEngineV2::StepKMC() {
  kengine_.StepKMC();
}

void UberEngineV2::DumpAll() {
  #ifdef DEBUG
  pengine_.Dump();
  fengine_.Dump();
  kengine_.Dump();
  #endif
}
