// Implementation for all pairs tracking scheme

#include "tracking_scheme_allpairs.h"

// Init funcitonality
void TrackingSchemeAllPairs::Init(int pModuleID,
                                  space_struct *pSpace,
                                  PotentialBase *pPotentialBase,
                                  std::vector<interaction_t> *pInteractions,
                                  std::vector<SpeciesBase*> *pSpecies,
                                  std::vector<Simple*> *pSimples,
                                  std::unordered_map<int, int> *pOIDMap,
                                  YAML::Node *pNode) {
  TrackingScheme::Init(pModuleID, pSpace, pPotentialBase, pInteractions, pSpecies, pSimples, pOIDMap, pNode);

  CreateTrackingScheme();
}

// Print functionality
void TrackingSchemeAllPairs::Print() {
  TrackingScheme::Print();
}

// Print statistics
void TrackingSchemeAllPairs::PrintStatistics() {
  GenerateStatistics();
  std::cout << "********\n";
  std::cout << name_ << std::endl;
  std::cout << "   {" << SIDToString(sid0_) << ", " << SIDToString(sid1_) << "}\n";
  std::cout << "   type: " << PtypeToString(type_) << std::endl;
  std::cout << "   nupdates: " << nupdates_ << std::endl;
  std::cout << "   avg time between updates: " << std::setprecision(8) << avg_update_time_/nupdates_ << " microseconds\n";
  std::cout << "   avg occpancy:             " << std::setprecision(8) << avg_occupancy_/nupdates_ << " particles\n";
}

// Create the all pairs tracking scheme
void TrackingSchemeAllPairs::CreateTrackingScheme() {
  GenerateAllPairs();
}

// Generate the interactions
void TrackingSchemeAllPairs::GenerateInteractions(bool pForceUpdate) {
  // Check to see if we have to update because the particle numbers changed, or something like that
  if (pForceUpdate) {
    GenerateStatistics(); // Generate before we update....
    LoadSimples();
    GenerateAllPairs();
  }
  interactions_->insert(interactions_->end(), m_interactions_.begin(), m_interactions_.end());
}

// Generate the all pairs stuff
// XXX FIXME maybe we don't need this, and can just loop over rigids
// and not do the checking at all?
void TrackingSchemeAllPairs::GenerateAllPairs() {
  // Clear the interactions
  nupdates_++;
  m_interactions_.clear();

  // Find the unique rigids
  unique_rids_->clear();
  for (int i = 0; i < nmsimples_; ++i) {
    unique_rids_->insert(m_simples_[i]->GetRID());
    mneighbors_[i].clear();
  }
  maxrigid_ = *(unique_rids_->rbegin());

  // clear out the kmc neighbor list
  if (type_ == ptype::kmc) {
    for (int i = 0; i < nsimples_; ++i) {
      neighbors_[i].clear();
    }
  }

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
    for (int idx = 0; idx < nmsimples_ - 1; ++idx) {
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

      for (int jdx = idx; jdx < nmsimples_; ++jdx) {
        if (idx == jdx) continue;
        auto p2 = m_simples_[jdx];
        int rid2 = p2->GetRID();
        if (rid1 == rid2) continue;

        // XXX FIXME this might not be optimal, but check anyway....
        auto sid0 = p1->GetSID();
        auto sid1 = p2->GetSID();
        if (!(sid0 == sid0_ && sid1 == sid1_) &&
            !(sid1 == sid0_ && sid0 == sid1_)) continue;

        // We are guranteed for an interaction
        if (rid_check_local->count(rid2)) {
          continue;
        }
        rid_check_local->insert(rid2);

        //// Create the interaction
        //interaction_t new_interaction;
        //new_interaction.idx_ = (*oid_position_map_)[p1->GetOID()];
        //new_interaction.jdx_ = (*oid_position_map_)[p2->GetOID()];
        //new_interaction.type_ = type_;
        //new_interaction.pot_ = pbase_;
        //new_interaction.kmc_track_module_ = moduleid_;
        //
        //// KMC specifics
        //if (type_ == ptype::kmc) {
        //  new_interaction.kmc_target_ = kmc_target_;
        //}
        //
        //#ifdef ENABLE_OPENMP
        //#pragma omp critical
        //#endif
        //{
        //  m_interactions_.push_back(new_interaction);
        //}
        neighbor_kmc_t new_neighbor;
        new_neighbor.idx_ = jdx;
        mneighbors_[idx].push_back(new_neighbor);

        // If necessary, build the main KMC neighbor list
        if (type_ == ptype::kmc) {
          neighbor_kmc_t new_neighbor_master;
          int nidx = (*oid_position_map_)[p1->GetOID()];
          int njdx = (*oid_position_map_)[p2->GetOID()];
          new_neighbor_master.idx_ = njdx; 
          neighbors_[nidx].push_back(new_neighbor_master);
        }
      } // for loop over second particle
    } // for loop over first particle
  } // omp parallel

  // Serialize the new interactions to the neighbor list
  for (int idx = 0; idx < nmsimples_; ++idx) {
    auto p1 = m_simples_[idx];
    nl_kmc_list *mlist = &mneighbors_[idx];
    for (auto nldx = mlist->begin(); nldx != mlist->end(); ++nldx) {
      auto p2 = m_simples_[nldx->idx_];
      interaction_t new_interaction;
      int nidx = (*oid_position_map_)[p1->GetOID()];
      int njdx = (*oid_position_map_)[p2->GetOID()];
      new_interaction.idx_ = nidx;
      new_interaction.jdx_ = njdx;
      new_interaction.type_ = type_;
      new_interaction.pot_ = pbase_;
      new_interaction.kmc_track_module_ = moduleid_;

      // KMC specifics
      if (type_ == ptype::kmc) {
        new_interaction.kmc_target_ = kmc_target_;

        // Create a bigger neighbor list so that we dno't have to
        // calculate it from the interactions later, just cache
        // it to the interaction
        neighbor_kmc_t *target_neighbor;
        for (auto nldx = neighbors_[nidx].begin(); nldx != neighbors_[nidx].end(); ++nldx) {
          if (nldx->idx_ == njdx) {
            target_neighbor = &(*nldx);
            break;
          }
        }
        new_interaction.neighbor_ = target_neighbor;
        if (debug_trace) {
          std::cout << "[KMC] Assigning [" << nidx << ", " << njdx << "] ->\n"
            << "   neighbor: " << new_interaction.neighbor_ << std::endl
            << "   original: " << &(neighbors_[nidx].back()) << std::endl
            << "   neighbor jdx: " << new_interaction.neighbor_->idx_ << std::endl;
        }
      }

      m_interactions_.push_back(new_interaction);
    }
  } //serialized add
}
