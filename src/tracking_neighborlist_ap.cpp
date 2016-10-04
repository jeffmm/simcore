// Implementation for neighbor list all pairs tracking

#include "tracking_neighborlist_ap.h"

void TrackingNeighborListAP::CreateSubstructure(double pRcut, nl_list** pNeighbors) {
  // Is there really anything to be done here?
  rcut_ = pRcut;
  rcut2_ = rcut_*rcut_;
  neighbors_ = (*pNeighbors);

  // Make sure we have a valid skin
  if (skin_ <= 0.0) {
    skin_ = 1.0;
  }

  skin2_ = skin_*skin_;
  rcs2_ = skin2_ + rcut2_ + 2*rcut_*skin_;
  half_skin2_ = 0.25*skin2_;
}

void TrackingNeighborListAP::UpdateRcut(double pRcut) {
  rcut_ = pRcut;
  rcut2_ = rcut_*rcut_;
  rcs2_ = skin2_ + rcut2_ + 2*rcut_*skin_;
}

void TrackingNeighborListAP::UpdateTracking(bool pForceUpdate) {
  nl_update_ = pForceUpdate;
  // Check for an update
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if (nl_update_)
      break;
    double dr2 = (*spec)->GetDrMax();
    nl_update_ = dr2 > half_skin2_;
  }
  if (nl_update_) {
    UpdateNeighborList();
  }
}

void TrackingNeighborListAP::UpdateNeighborList() {
  if (debug_trace) {
    std::cout << "Updating NeighborListAP\n";
  }
  nupdates_++;
  unique_rids_->clear();
  for (int i = 0; i < nsimples_; ++i) {
    neighbors_[i].clear();
    unique_rids_->insert((*simples_)[i]->GetRID());
  }
  maxrigid_ = *(unique_rids_->rbegin());

  // Call the main update routine
  AllPairsUpdate2();

  // Reset the accumulators
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->ZeroDr();
  }
}

/*
void TrackingNeighborListAP::UpdateNeighborList() {
  if (debug_trace) {
    std::cout << "Updating NeighborListAP\n";
  }
  nupdates_++;
  for (int i = 0; i < nsimples_; ++i) {
    neighbors_[i].clear();
  }

  // Call the main update routine
  AllPairsUpdate();

  // Reset the accumulators
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->ZeroDr();
  }
}
*/

void TrackingNeighborListAP::AllPairsUpdate2() {
  // Loop over all pairs to build the neighbor list, O(N^2)

  // Clear the rid self check
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
    for (int idx = 0; idx < nsimples_; ++idx) {
      auto p1 = (*simples_)[idx];
      int rid1 = p1->GetRID();

      bool should_exit = false;
      #ifdef ENABLE_OPENMP
      #pragma omp critical
      #endif
      {
        should_exit = (*rid_self_check_)[rid1];
        (*rid_self_check_)[rid1] = true;
      }
      if (should_exit) continue;

      // Get the threadlocal rid_check_local_
      std::unordered_set<int>* rid_check_local = rid_check_local_[tid];
      rid_check_local->clear();

      for (int jdx = 0; jdx < nsimples_; ++jdx) {
        if (idx == jdx) continue;
        auto p2 = (*simples_)[jdx];
        int rid2 = p2->GetRID();
        if (rid1 == rid2) continue;

        // Check if we've already seen this rid for ourselves
        if (rid_check_local->count(rid2)) {
          continue;
        }

        // Check if there is even an interaction
        PotentialBase *pot1 = potentials_->GetPotentialExternal(p1->GetSID(), p2->GetSID());
        PotentialBase *pot2 = potentials_->GetPotentialInternal(p1->GetOID(), p2->GetOID());
        if ((pot1 == nullptr) && (pot2 == nullptr)) continue;

        // Minimum distance
        interactionmindist idm;
        MinimumDistance(p1, p2, idm, ndim_, nperiodic_, space_);

        if (idm.dr_mag2 < rcs2_) {
          neighbor_t new_neighbor;
          new_neighbor.idx_ = jdx;
          new_neighbor.rid_me_ = rid1;
          new_neighbor.rid_you_ = rid2;
          neighbors_[idx].push_back(new_neighbor);
          rid_check_local->insert(rid2);
        }
      }
    } // pragma omp for

    #ifdef ENABLE_OPENMP
    #pragma omp barrier
    #endif
  } // omp parallel
}

void TrackingNeighborListAP::AllPairsUpdate() {
  // Loop over all pairs to build the neighbor list, O(N^2)
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
    //std::set<std::pair<int, int>> rid_interactions;
    rid_interactions_[tid].clear();
    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int idx = 0; idx < nsimples_; ++idx) {
      for (int jdx = 0; jdx < nsimples_; ++jdx) {
        if (idx == jdx) continue;
        auto p1 = (*simples_)[idx];
        auto p2 = (*simples_)[jdx];
        if (p1->GetOID() == 0) dummy_function();
        int rid1 = p1->GetRID();
        int rid2 = p2->GetRID();
        // Only insert object if rigid ID is unique
        // std::set only allows insertion of unique elements
        // Check if insertion fails because RID is not unique
        //if (!rid_interactions.insert(std::make_pair(rid1, rid2)).second)
        //  continue;
        if (!rid_interactions_[tid].insert(std::make_pair(rid1, rid2)).second)
          continue;
        // XXX Don't exclude self interactions of composite objects
        // needed for filaments and xlinks internal properties, things like
        // that
        if (p1->GetRID() == p2->GetRID()) continue;

        // Minimum distance
        interactionmindist idm;
        MinimumDistance(p1, p2, idm, ndim_, nperiodic_, space_);

        if (idm.dr_mag2 < rcs2_) {
          neighbor_t new_neighbor;
          new_neighbor.idx_ = jdx;
          new_neighbor.rid_me_ = rid1;
          new_neighbor.rid_you_ = rid2;
          neighbors_[idx].push_back(new_neighbor);
        }
      }
    } // pragma omp for
  } // pragma omp parallel
  #ifdef ENABLE_OPENMP
  {
    // Clean up any existing rid-rid interactions
    //std::set< std::pair<int, int> > rid_interactions;
    rid_interactions_[nthreads_].clear();
    for (int idx=0; idx<nsimples_; ++idx) {
      for (auto nb = neighbors_[idx].begin(); nb != neighbors_[idx].end();) {
        //if (!rid_interactions.insert(std::make_pair(nb->rid_me_,nb->rid_you_)).second) {
        if (!rid_interactions_[nthreads_].insert(std::make_pair(nb->rid_me_, nb->rid_you_)).second) {
            //if (debug_trace)
              //printf("Removing interaction pair (%d, %d)\n",nb->rid_you_,nb->rid_me_);
          neighbors_[idx].erase(nb);
        }
        else {
          ++nb;
        }
      }
    }

  }
  #endif
}

void TrackingNeighborListAP::print() {
  printf("\t{rcut: %2.2f}, {skin: %2.2f} = {rcs2: %2.2f}, {half_skin2: %2.2f}\n", rcut_, skin_, rcs2_, half_skin2_);
}

void TrackingNeighborListAP::dump() {
  printf("\tneighbor list-all pairs has nothing interesting to dump\n");
}
