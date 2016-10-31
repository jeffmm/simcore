// Implementation of tracking scheme base

#include "tracking_scheme.h"

// Initialize
void TrackingScheme::Init(int pModuleID,
                          space_struct *pSpace,
                          PotentialBase *pPotentialBase,
                          std::vector<interaction_t> *pInteractions,
                          std::vector<SpeciesBase*> *pSpecies,
                          std::vector<Simple*> *pSimples,
                          std::unordered_map<int, int> *pOIDMap,
                          YAML::Node *pNode) {

  space_ = pSpace;
  pbase_ = pPotentialBase;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  interactions_ = pInteractions;
  simples_ = pSimples;
  oid_position_map_ = pOIDMap;
  moduleid_ = pModuleID;

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

  YAML::Node node = *pNode;

  // SIDs for this interaction
  std::string sid0s = node["sid1"].as<std::string>();
  std::string sid1s = node["sid2"].as<std::string>();
  sid0_ = StringToSID(sid0s);
  sid1_ = StringToSID(sid1s);

  symmetric_ = (sid0_ == sid1_) ? true : false;

  // Type of interaction
  std::string types = node["type"].as<std::string>();
  type_ = StringToPtype(types);

  // KMC target if there is one
  if (node["kmc_target"]) {
    std::string kmc_targets = node["kmc_target"].as<std::string>();
    kmc_target_ = StringToSID(kmc_targets);
  }

  // Locate the two species we're interested in
  for (auto spec = pSpecies->begin(); spec != pSpecies->end(); ++spec) {
    if ((*spec)->GetSID() == sid0_) {
      spec0_ = (*spec);
    }
    if ((*spec)->GetSID() == sid1_) {
      spec1_ = (*spec);
    }
  }

  if (spec0_ == nullptr) {
    std::cout << "Coudln't find first species: " << sid0s << ", exiting\n";
    exit(1);
  }
  if (spec1_ == nullptr) {
    std::cout << "Coudln't find second species: " << sid1s << ", exiting\n";
    exit(1);
  }

  LoadSimples();

  // Generate the rigid checkers
  rid_self_check_ = new std::vector<bool>();
  unique_rids_ = new std::set<int>();
  rid_check_local_ = new std::unordered_set<int>*[nthreads_];
  for (int ithread = 0; ithread < nthreads_; ++ithread) {
    rid_check_local_[ithread] = new std::unordered_set<int>();
  }
}

// Print out the information for this tracking scheme
void TrackingScheme::Print() {
  std::cout << name_ << " Tracking Scheme [" << moduleid_ << "]\n";
  std::cout << "   symmetric: " << (symmetric_ ? "true" : "false") << std::endl;
  std::cout << "   potential: ->\n";
  std::cout << "   {" << SIDToString(sid0_) << ", " << SIDToString(sid1_) << "} : ";
  pbase_->Print();
  std::cout << "   type: " << PtypeToString(type_) << std::endl;
  std::cout << "   my nsimples: " << nmsimples_ << std::endl;
  if (type_ == ptype::kmc) {
    std::cout << "   kmc_target: " << SIDToString(kmc_target_) << std::endl;
  }
}

// Dump the internal neighbor list stuff
void TrackingScheme::Dump() {
  #ifdef DEBUG
  std::cout << name_ << " Tracking Scheme [" << moduleid_ << "] -> Dump\n";
  std::cout << "--------\n";
  for (int idx = 0; idx < nmsimples_; ++idx) {
    std::cout << "[" << idx << "] -> [";
    nl_kmc_list *mlist = &mneighbors_[idx];
    for (auto nldx = mlist->begin(); nldx != mlist->end(); ++nldx) {
      std::cout << nldx->idx_ << ", ";
    }
    std::cout << "]\n";
  }
  std::cout << "--------\n";
  #endif
}

// Load the simples symmetric or asymmetric
void TrackingScheme::LoadSimples() {
  if (debug_trace) {
    std::cout << "TrackingScheme::LoadSimples\n";
  }
  nsimples_ = (int)simples_->size();

  std::vector<Simple*> sim_vec0 = spec0_->GetSimples();
  std::vector<Simple*> sim_vec1 = spec1_->GetSimples();

  m_simples_.clear();
  nmsimples0_ = (int)sim_vec0.size();
  nmsimples1_ = (int)sim_vec1.size();
  m_simples_.insert(m_simples_.end(), sim_vec0.begin(), sim_vec0.end());
  if (!symmetric_) {
    m_simples_.insert(m_simples_.end(), sim_vec1.begin(), sim_vec1.end());
  }
  nmsimples_ = (int)m_simples_.size();

  // Set up the neighbor list
  if (mneighbors_ != nullptr) {
    delete[] mneighbors_;
  }
  mneighbors_ = new nl_kmc_list[nmsimples_];

  // Ugh, figure out the OID stuff
  m_oid_position_map_.clear();
  for (int i = 0; i < nmsimples_; ++i) {
    auto part = m_simples_[i];
    int oid = part->GetOID();
    m_oid_position_map_[oid] = i;
  }
}

// Generate the statistics
void TrackingScheme::GenerateStatistics() {
  // Generate avg occupancy of neighbor list, etc
  // Get the time
  this_time_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> elapsed_microseconds = this_time_ - last_time_;
  if (debug_trace) {
    std::cout << "[TrackingScheme] Update Elapsed: " << std::setprecision(8) << elapsed_microseconds.count() << std::endl;
  }
  last_time_ = this_time_;
  avg_update_time_ += elapsed_microseconds.count();

  // Get the occupancy
  int nneighbors = 0;
  for (int idx = 0; idx < nmsimples_; ++idx) {
    nl_kmc_list *mlist = &mneighbors_[idx];
    nneighbors += (int)mlist->size();
  }
  // XXX FIXME switch out for the nsteps that elapsed, and do over all steps of the system
  avg_occupancy_ += (double)nneighbors/(double)nmsimples_;
  if (debug_trace) {
    std::cout << "[TrackingScheme] Avg Occupancy: " << std::setprecision(8) << (double)nneighbors/(double)nmsimples_ << std::endl;
  }
}
