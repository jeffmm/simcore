#include "particle_engine.h"

#include "tracking_scheme_allpairs.h"
#include "tracking_scheme_nlap.h"

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
  if (debug_trace) {
    std::cout << "ParticleEngine::Init\n";
  }
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

}

void ParticleEngine::Print() {
  std::cout << "********\n";
  std::cout << "Particle Tracking Engine -> Schemes\n";
  for (int i = 0; i < (int)tracking_.size(); ++i) {
    std::cout << "[" << i << "] : ";
    tracking_[i]->Print();
  }
  PrintPotentials();
}

void ParticleEngine::PrintPotentials() {
  std::cout << "********\n";
  std::cout << "Particle Tracking Engine -> Potentials\n";
  potentials_.Print();
}

void ParticleEngine::Dump() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "ParticleEngine::Dump\n";
  //DumpSimples();
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
  REGISTER_SCHEME(TrackingSchemeNeighborListAllPairs,neighborlistallpairs);
}

// Load our simples array
void ParticleEngine::LoadSimples() {
  if (debug_trace) {
    std::cout << "ParticleEngine::LoadSimples\n";
  }
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
  if (debug_trace) {
    std::cout << "ParticleEngine::CreateTracking\n";
  }
  // Load the potential file
  char *fname = params_->potfile;
  node_ = YAML::LoadFile(fname);
  npots_ = (int)node_["potentials"].size();
  for (int ipot = 0; ipot < npots_; ++ipot) {
    // Add the potential to the potential manager, and get the entry
    YAML::Node subnode = node_["potentials"][ipot];
    int potidx = potentials_.AddPotential(&subnode);

    std::string types   = node_["potentials"][ipot]["type"].as<std::string>();
    // Determine what kind of potential we have, as it affects if we attach to something or not
    if (types.compare("external") == 0) {
      CreateExternalPotential(&subnode, potidx);
    } else if (types.compare("internal") == 0) {
      CreateInternalPotential(&subnode, potidx); 
    } else if (types.compare("boundary") == 0) {
      CreateBoundaryPotential(&subnode, potidx);
    } else if (types.compare("tether") == 0) {
      CreateTetherPotential(&subnode, potidx);
    } else {
      std::cout << "Please specify a potential type we understand, exiting\n";
      exit(1);
    }
  }
}

// Create an external potential
void ParticleEngine::CreateExternalPotential(YAML::Node *subnode, int potidx) {
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
}

// Create the internal potential
void ParticleEngine::CreateInternalPotential(YAML::Node *subnode, int potidx) {
  YAML::Node node = *subnode;
  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // There is no scheme, just raw information
  std::string sids = node["sid"].as<std::string>();
  SID sid = StringToSID(sids);
  SpeciesBase *sit;
  for (auto msit = species_->begin(); msit != species_->end(); ++msit) {
    if ((*msit)->GetSID() == sid) {
      sit = (*msit);
      break;
    }
  }
  auto internal_pairs = sit->GetInternalPairs();
  for (auto ipair = internal_pairs.begin(); ipair != internal_pairs.end(); ++ipair) {
    interaction_t new_interaction;
    new_interaction.idx_ = oid_position_map_[ipair->first];
    new_interaction.jdx_ = oid_position_map_[ipair->second];
    new_interaction.type_ = ptype::internal;
    new_interaction.pot_ = mypot;
    m_internal_interactions_.push_back(new_interaction);
  }
}

// Creat boundary potentials
void ParticleEngine::CreateBoundaryPotential(YAML::Node *subnode, int potidx) {
  YAML::Node node = *subnode;
  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // No scheme, just attach to boundary stuff
  std::string sids = node["sid"].as<std::string>();
  SID sid = StringToSID(sids);
  SpeciesBase *sit;
  for (auto msit = species_->begin(); msit != species_->end(); ++msit) {
    if ((*msit)->GetSID() == sid) {
      sit = (*msit);
      break;
    }
  }

  auto spec_simples = sit->GetSimples();
  // Loop over the spec_members, getting unique element 0 simples from each
  // member
  // Each RID can interact with a boundary
  int nsimp = (int)spec_simples.size();
  int last_rid = -1;
  for (int idx = 0; idx < nsimp; ++idx) {
    auto part = spec_simples[idx];
    if (part->GetRID() != last_rid) {
      last_rid = part->GetRID();
      interaction_t new_interaction;
      new_interaction.idx_ = oid_position_map_[part->GetOID()];
      new_interaction.jdx_ = -1;
      new_interaction.type_ = ptype::boundary;
      new_interaction.pot_ = mypot;
      m_boundary_interactions_.push_back(new_interaction);
    } else {
      // Just continue, not a unique RID for this
    }
  }
}

// Create a tether potential
void ParticleEngine::CreateTetherPotential(YAML::Node *subnode, int potidx) {
  YAML::Node node = *subnode;
  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // No scheme, but have to use the anchor list to attach
  for (auto ait = anchors_->begin(); ait != anchors_->end(); ++ait) {
    // What anchor are we on
    std::cout << "Anchor: " << ait->first << std::endl;
    std::vector<anchor_t>* avec = &(ait->second);
    for (auto ait2 = avec->begin(); ait2 != avec->end(); ++ait2) {
      // Get the pointer to the anchor
      anchor_t *manchor = &(*ait2);
      std::cout << "   -> {" << manchor->idx_base_ << ", " << manchor->idx_other_ << "}\n";
      interaction_t new_interaction;
      new_interaction.idx_ = oid_position_map_[manchor->idx_base_];
      new_interaction.jdx_ = oid_position_map_[manchor->idx_other_];
      new_interaction.type_ = ptype::tether;
      new_interaction.pot_ = mypot;
      new_interaction.anchor_ = manchor;
      m_tether_interactions_.push_back(new_interaction);
    }
  }
}

// Create a special KMC tracking node
TrackingScheme* ParticleEngine::CreateKMCExternal(YAML::Node *subnode) {
  YAML::Node node = *subnode;
  std::string types = node["type"].as<std::string>();
  if (types.compare("kmc") != 0) {
    std::cout << "Attempting a KMC load potential with " << types << std::endl;
    exit(1);
  }

  // Create the potential(s)
  int potidx = potentials_.AddPotential(&node);

  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // Get the listed scheme
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

  return scheme;
}

// Create a special KMC internal set of interactions (eg xlinks)
PotentialBase *ParticleEngine::CreateKMCInternal(YAML::Node *subnode) {
  YAML::Node node = *subnode;
  std::string types = node["type"].as<std::string>();
  if (types.compare("internal") != 0) {
    std::cout << "Attempting a KMC internal load potential with " << types << std::endl;
    exit(1);
  }

  // Create the potential
  int potidx = potentials_.AddPotential(&node);

  PotentialBase *mypot = potentials_.GetPotential(potidx);

  // Internal, so there is no scheme
  CreateInternalPotential(&node, potidx);

  return mypot;
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
void ParticleEngine::UpdateInteractions(bool pForceUpdate) {
  // Clear interactions, will reload with cached version if no update needed, 
  // but this is handled at the scheme level
  interactions_->clear();
  // Check for an update that requires rebuilding the simples, etc
  trigger_update_ = pForceUpdate;
  CheckTriggerUpdate();
  // External and kmc interactions
  for (auto ixs = tracking_.begin(); ixs != tracking_.end(); ++ixs) {
    (*ixs)->GenerateInteractions(trigger_update_);
  }
  // Add the internal interactions
  interactions_->insert(interactions_->end(), m_internal_interactions_.begin(), m_internal_interactions_.end());
  // Boundary interactions
  interactions_->insert(interactions_->end(), m_boundary_interactions_.begin(), m_boundary_interactions_.end());
  // Tethers
  interactions_->insert(interactions_->end(), m_tether_interactions_.begin(), m_tether_interactions_.end());
}

