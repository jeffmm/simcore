// Implementation for all pairs tracking scheme

#include "tracking_scheme_allpairs.h"

// Init funcitonality
void TrackingSchemeAllPairs::Init(space_struct *pSpace,
                                  PotentialBase *pPotentialBase,
                                  std::vector<interaction_t> *pInteractions,
                                  std::vector<SpeciesBase*> *pSpecies,
                                  std::vector<Simple*> *pSimples,
                                  std::unordered_map<int, int> *pOIDMap,
                                  YAML::Node *pNode) {
  TrackingScheme::Init(pSpace, pPotentialBase, pInteractions, pSpecies, pSimples, pOIDMap, pNode);

  CreateTrackingScheme();
}

// Print functionality
void TrackingSchemeAllPairs::Print() {
  TrackingScheme::Print();
}

// Create the all pairs tracking scheme
void TrackingSchemeAllPairs::CreateTrackingScheme() {
  std::cout << "TrackingSchemeAllPairs CreateTrackingScheme\n";
  GenerateAllPairs();
}

// Generate the interactions
void TrackingSchemeAllPairs::GenerateInteractions(bool pForceUpdate) {
  std::cout << "TrackingSchemeAllPairs GenerateInteractions\n";
  // Check to see if we have to update because the particle numbers changed, or something like that
  //CheckTriggerUpdate();
  if (pForceUpdate) {
    LoadSimples();
    GenerateAllPairs();
  }
  interactions_->insert(interactions_->end(), m_interactions_.begin(), m_interactions_.end());
}

// Generate the all pairs stuff
void TrackingSchemeAllPairs::GenerateAllPairs() {
  std::cout << "TrackingSchemeAllPairs GenerateAllPairs\n";
  // Clear the interactions
  nupdates_++;
  m_interactions_.clear();

  // Find the unique rigids
  unique_rids_->clear();
  for (int i = 0; i < nmsimples_; ++i) {
    unique_rids_->insert(m_simples_[i]->GetRID());
  }
  maxrigid_ = *(unique_rids_->rbegin());

  // Loop over particles, generate the interactions
  rid_self_check_->clear();
  rid_self_check_->resize(maxrigid_+1);
  std::fill(rid_self_check_->begin(), rid_self_check_->end(), false);

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    int tid = 0;
    #ifdef ENABLE_OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int idx = 0; idx < nmsimples_; ++idx) {
      auto p1 = m_simples_[idx];
      int rid1 = p1->GetRID();

      bool should_exit = false;
      #ifdef ENABLE_OPENMP
      #pragma omp critical
      #endif
      {
        should_exit = (*rid_self_check_)[rid1];
        (*rid_self_check_)[rid1] = true;
      }
      if (should_exit) {
        continue;
      }

      std::unordered_set<int>* rid_check_local = rid_check_local_[tid];
      rid_check_local->clear();

      for (int jdx = 0; jdx < nmsimples_; ++jdx) {
        if (idx == jdx) continue;
        auto p2 = m_simples_[jdx];
        int rid2 = p2->GetRID();
        if (rid1 == rid2) continue;

        // We are guranteed for an interaction
        if (rid_check_local->count(rid2)) {
          continue;
        }
        rid_check_local->insert(rid2);

        // Create the interaction
        interaction_t new_interaction;
        new_interaction.idx_ = (*oid_position_map_)[p1->GetOID()];
        new_interaction.jdx_ = (*oid_position_map_)[p2->GetOID()];
        new_interaction.type_ = type_;
        new_interaction.pot = pbase_;
        
        #ifdef ENABLE_OPENMP
        #pragma omp critical
        #endif
        {
          m_interactions_.push_back(new_interaction);
        }



      } // for loop over second particle
    } // for loop over first particle
  } // omp parallel
}
