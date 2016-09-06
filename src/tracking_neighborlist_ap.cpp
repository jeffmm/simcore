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

void TrackingNeighborListAP::AllPairsUpdate() {
  // Loop over all pairs to build the neighbor list, O(N^2)
  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    std::set<std::pair<int, int>> rid_interactions;
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
        if (!rid_interactions.insert(std::make_pair(rid1, rid2)).second)
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
    std::set< std::pair<int, int> > rid_interactions;
    for (int idx=0; idx<nsimples_; ++idx) {
      for (auto nb = neighbors_[idx].begin(); nb != neighbors_[idx].end();) {
        if (!rid_interactions.insert(std::make_pair(nb->rid_me_,nb->rid_you_)).second) {
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
