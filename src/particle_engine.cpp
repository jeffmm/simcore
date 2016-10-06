#include "particle_engine.h"

#include "tracking_scheme_allpairs.h"

#define REGISTER_SCHEME(n,m) scheme_factory_.register_class<n>(#m);

// Destructor
ParticleEngine::~ParticleEngine() {
  std::cout << "********\n";
  std::cout << "Particle Engine - Final Statistics\n";
  for (auto ixs = tracking_.begin(); ixs != tracking_.end(); ++ixs) {
    (*ixs)->PrintStatistics();
    delete (*ixs);
  }
}

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
  potentials_.Init(space_, anchors_);

  RegisterSchemes();

  LoadSimples();

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

void ParticleEngine::Dump() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "ParticleEngine::Dump\n";
  DumpSimples();
  //DumpInteractions();
  #endif
}

void ParticleEngine::DumpInteractions() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "DumpInteractions\n";
  for (int ixs = 0; ixs < interactions_->size(); ++ixs) {
    auto mixs = (*interactions_)[ixs];
    std::cout << "[" << ixs << "] {" << mixs.idx_ << " -> " << mixs.jdx_ << "}, type: "
      << PtypeToString(mixs.type_) << std::endl;
  }
  #endif
}

void ParticleEngine::DumpSimples() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "DumpSimples\n";
  for (int idx = 0; idx < (int)simples_.size(); ++idx) {
    auto part = simples_[idx];
    part->Dump();
  }
  #endif
}

// Register the tracking schemes we know about
void ParticleEngine::RegisterSchemes() {
  REGISTER_SCHEME(TrackingSchemeAllPairs,allpairs);
}

// Load our simples array
void ParticleEngine::LoadSimples() {
  nsys_ = (int)species_->size();
  simples_.clear();
  for (int ispec = 0; ispec < nsys_; ++ispec) {
    std::vector<Simple*> sim_vec = (*species_)[ispec]->GetSimples();
    simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
  }
  nsimples_ = (int)simples_.size();

  // Ugh, figure out the OID stuff
  oid_position_map_.clear();
  for (int i = 0; i < nsimples_; ++i) {
    auto part = simples_[i];
    int oid = part->GetOID();
    oid_position_map_[oid] = i;
  }
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
  scheme->Init((int)tracking_.size(),
               space_,
               mypot,
               interactions_,
               species_,
               &simples_,
               &oid_position_map_,
               &node);
  tracking_.push_back(scheme);

  std::cout << "Particle Engine CreateExternalPotential done\n";
}

// Check if a global update has been triggered
void ParticleEngine::CheckTriggerUpdate() {
  int simples_count = 0;
  for (auto ispec = species_->begin(); ispec != species_->end(); ++ispec) {
    simples_count += (*ispec)->GetCount();
  }
  if (simples_count != nsimples_) {
    if (debug_trace) {
      std::cout << "ParticleEngine nsimples changed, forcing update\n";
    }
    trigger_update_ = true;
  }
  for (int ispec = 0; ispec < nsys_; ++ispec) {
    if ((*species_)[ispec]->GetUpdate()) {
      if (debug_trace) {
        std::cout << "ParticleEngine species should update, forcing update\n";
      }
      (*species_)[ispec]->SetUpdate(false);
      trigger_update_ = true;
    }
  }

  if (trigger_update_) {
    LoadSimples();
  }
}

// Generate the interactions
void ParticleEngine::UpdateInteractions() {
  // Clear interactions, will reload with cached version if no update needed, 
  // but this is handled at the scheme level
  interactions_->clear();
  // Check for an update that requires rebuilding the simples, etc
  trigger_update_ = false;
  CheckTriggerUpdate();
  // External and kmc interactions
  for (auto ixs = tracking_.begin(); ixs != tracking_.end(); ++ixs) {
    (*ixs)->GenerateInteractions(trigger_update_);
  }

}

