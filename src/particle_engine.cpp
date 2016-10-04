#include "particle_engine.h"

#include "tracking_scheme_allpairs.h"

#define REGISTER_SCHEME(n,m) scheme_factory_.register_class<n>(#m);

// Pass in the main system properties information
void ParticleEngine::Init(system_parameters *pParams,
                          space_struct *pSpace,
                          std::vector<SpeciesBase*> *pSpecies,
                          al_set *pAnchors,
                          std::vector<interaction_t> *pInteractions,
                          long seed) {
  std::cout << "Particle Engine Init\n";
  params_ = pParams;
  space_ = pSpace;
  species_ = pSpecies;
  anchors_ = pAnchors;
  interactions_ = pInteractions;
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

  // Init the potential manager
  potentials_.Init(species_, space_, anchors_);

  RegisterSchemes();

  std::cout << "Particle Engine Init done\n";
}

void ParticleEngine::Print() {
  std::cout << "********\n";
  std::cout << "Particle Tracking Engine ->\n";
  for (int i = 0; i < (int)tracking_.size(); ++i) {
    std::cout << "[" << i << "] : ";
    tracking_[i]->Print();
  }
}

// Register the tracking schemes we know about
void ParticleEngine::RegisterSchemes() {
  REGISTER_SCHEME(TrackingSchemeAllPairs,allpairs);
}

// Create all of the tracking information
void ParticleEngine::CreateTracking() {
  std::cout << "Particle Engine CreateTracking\n";

  // Load the potential file
  char *fname = params_->potfile;
  std::cout << "Loading potentials from " << fname << std::endl;
  node_ = YAML::LoadFile(fname);
  npots_ = (int)node_["potentials"].size();
  std::cout << "Found " << npots_ << " potentials\n";
  for (int ipot = 0; ipot < npots_; ++ipot) {
    // Add the potential to the potential manager, and get the entry
    YAML::Node subnode = node_["potentials"][ipot];
    int potidx = potentials_.AddPotential(&subnode);
    std::cout << "Put pot in idx: " << potidx << std::endl;

    std::string types   = node_["potentials"][ipot]["type"].as<std::string>();
    // Determine what kind of potential we have, as it affects if we attach to something or not
    if (types.compare("external") == 0 ||
        types.compare("kmc") == 0) {
      std::cout << "External or kmc, need tracking\n";
      CreateExternalPotential(&subnode, potidx);
    } else if (types.compare("internal") == 0 ||
               types.compare("boundary") == 0 ||
               types.compare("tether") == 0) {
      std::cout << "Found one of the other types\n";
    } else {
      std::cout << "Please specify a potential type we understand, exiting\n";
      exit(1);
    }
  }

  //potentials_.Print();

  std::cout << "Particle Engine CreateTracking done\n";
}

void ParticleEngine::CreateExternalPotential(YAML::Node *subnode, int potidx) {
  std::cout << "Particle Engine CreateExternalPotential\n";

  YAML::Node node = *subnode;
  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // Look for the scheme type listed
  std::string schemes = node["scheme"].as<std::string>();
  TrackingScheme *scheme = (TrackingScheme*) scheme_factory_.construct(schemes);

  if (!scheme) {
    std::cout << "Scheme " << schemes << " not found, exiting\n";
    exit(1);
  }
  scheme->Init(space_, mypot, interactions_, &node);


  //scheme->Print();
  tracking_.push_back(scheme);


  std::cout << "Particle Engine CreateExternalPotential done\n";
}
