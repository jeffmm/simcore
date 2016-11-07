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

  // Set the engine to calculate neighbor lists for kmc particles on the first go around
  pengine_.UpdateInteractions(true);
  pengine_.SetTriggerUpdate(true);

  fengine_.Interact();
  kengine_.PreGenerateNeighbors();

  pengine_.Print();
  kengine_.Print();
  pengine_.Dump();
  fengine_.Dump();
  kengine_.Dump();

  last_time_ = std::chrono::high_resolution_clock::now();
  ndatapoints_ = 0;
  ninteractions_ = 0;
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

void UberEngineV2::GenerateStatistics(int istep) {
  // Generate the statistics for this step
  // Compute the elapsed time
  if (istep % 100 != 0) return;
  ndatapoints_++;
  this_time_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> elapsed_microseconds = this_time_ - last_time_;
  last_time_ = this_time_;

  // Compute interaction quantities
  ninteractions_ += (int)interactions_.size();

  if (debug_trace) {
    std::cout << "--------\n";
    std::cout << "Generate Statistics ->\n";
    std::cout << "Elapsed time: " << std::setprecision(8) << elapsed_microseconds.count() << " microseconds\n";
    std::cout << "N interactions (this delta): " << interactions_.size() << std::endl;
    std::cout << "N interactions avg:          " << std::setprecision(8) << (double)ninteractions_/(int)ndatapoints_ << std::endl;
  }

}

void UberEngineV2::PrintStatistics() {
  std::cout << "********\n";
  std::cout << "UberEngineV2 - Final Statistics\n";
  std::cout << "N interactions avg: " << std::setprecision(8) << (double)ninteractions_/(int)ndatapoints_ << std::endl;
}

void UberEngineV2::PrepOutputs() {
  kengine_.PrepOutputs();
}

void UberEngineV2::WriteOutputs(int istep) {
  kengine_.WriteOutputs(istep);
}
