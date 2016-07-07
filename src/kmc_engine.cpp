// Implementation for kinetic monte carlo engine

#include "kmc_engine.h"

void kmcEngine::Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, ParticleTracking *pTracking) {
  space_ = pSpace;
  species_ = pSpecies;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;

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
}

// Step the kmc engine forward one
void kmcEngine::StepKMC() {
  PrepKMC();

  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->StepKMC();
    }
  }
}

// Prepare and update probabilities of the kmc engine
void kmcEngine::PrepKMC() {
  // Just ask each species what is needs to do to prep
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->PrepKMC();
    }
  }
}
