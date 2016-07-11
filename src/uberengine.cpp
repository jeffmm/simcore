#include <chrono>

#include "UberEngine.h"
#include "object.h"

// Pass in the main system properties information
void UberEngine::Init(system_parameters *pParams, space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, long seed) {
  params_ = pParams;
  space_ = pSpace;
  species_ = pSpecies;
  n_dim_ = space_->n_dim;
  n_periodic_ = space_->n_periodic;
  max_overlap_ = params_->max_overlap;
  draw_flag_ = params_->draw_interactions;
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

  // Switch for interaction substructure type
  printf("********\n");
  printf("Creating UberEngine\n");
  switch (params_->ftype) {
    case 0:
        printf("Must specify a tracking substructure, exiting!\n");
        exit(1);
    case 1:
        force_type_ = FTYPE::allpairs;
        skin_ = 0.0;
        break;
    case 2:
        force_type_ = FTYPE::microcells;
        skin_ = 0.0;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    case 3:
        force_type_ = FTYPE::cells;
        skin_ = 0.0;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    case 4:
        force_type_ = FTYPE::neighborallpairs;
        skin_ = params_->masterskin;
        break;
    case 5:
        force_type_ = FTYPE::neighborcells;
        skin_ = params_->masterskin;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    default:
        printf("Must specify a force substructure, exiting!\n");
        break;
   }
  InitPotentials();

  // Initialize the particle tracking
  tracking_.Init(space_, species_, skin_, force_type_);
  tracking_.InitPotentials(&potentials_);
  tracking_.InitTracking();
  tracking_.CheckOverlaps(max_overlap_);
  tracking_.Print();

  // Initialize the interaction engine
  fengine_.Init(space_, &tracking_, skin_);
  fengine_.InitPotentials(&potentials_);
  fengine_.InitMP();
  fengine_.Print();

  // Initialize the kmc engine
  kengine_.Init(space_, species_, &tracking_, gsl_rng_get(rng_.r));
  kengine_.InitMP();

  // Run one step to make sure that we're all good
  tracking_.UpdateTracking(true);
  fengine_.Interact();

}

void UberEngine::InitPotentials() {
  // Ask the potential manager to parse the potentials file
  potentials_.Init(space_, params_->potfile);
}

void UberEngine::DumpAll() {
    // Dump all the particles and their positions, forces, energy (2d)
    #ifdef DEBUG
    tracking_.Dump();
    fengine_.Dump();
    kengine_.Dump();
    #endif
}

void UberEngine::InteractMP() {
  tracking_.UpdateTracking();
  fengine_.Interact();
}

void UberEngine::StepKMC() {
  kengine_.StepKMC();
}

void UberEngine::Draw(std::vector<graph_struct*> * graph_array) {
  std::cerr << "Drawing forces currently doesn't work, exiting\n";
  exit(1);
  if (!draw_)
    return;
  for (auto it=draw_array_.begin(); it!=draw_array_.end(); ++it)
    graph_array->push_back(&(*it));
}


