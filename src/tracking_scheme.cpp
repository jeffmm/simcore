// Implementation of tracking scheme base

#include "tracking_scheme.h"

// Initialize
void TrackingScheme::Init(space_struct *pSpace,
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

  // Type of interaction
  std::string types = node["type"].as<std::string>();
  type_ = StringToPtype(types);

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
  nsimples_ = (int)simples_->size();

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
  std::cout << name_ << " Tracking Scheme" << std::endl;
  std::cout << "   potential: ->\n";
  std::cout << "   {" << SIDToString(sid0_) << ", " << SIDToString(sid1_) << "} : ";
  pbase_->Print();
  std::cout << "   type: " << PtypeToString(type_) << std::endl;
  std::cout << "   my nsimples: " << nmsimples_ << std::endl;
}

// Overall load simples functionality
void TrackingScheme::LoadSimples() {
  std::cout << "TrackingScheme Load Simples\n";
  std::vector<Simple*> sim_vec0 = spec0_->GetSimples();
  std::vector<Simple*> sim_vec1 = spec1_->GetSimples();

  m_simples_.insert(m_simples_.end(), sim_vec0.begin(), sim_vec0.end());
  // Avoiding adding twice if we're the same species
  if (sid0_ != sid1_) {
    m_simples_.insert(m_simples_.end(), sim_vec1.begin(), sim_vec1.end());
  }

  nmsimples_ = (int)m_simples_.size();
}

