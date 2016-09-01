// Implementation for all pairs (brute) force tracking

#include "tracking_allpairs.h"

void TrackingAllPairs::CreateSubstructure(double pRcut, nl_list** pNeighbors) {
  // Is there really anything to be done here?
  rcut_ = pRcut;
  neighbors_ = (*pNeighbors);
}

void TrackingAllPairs::UpdateRcut(double pRcut) {
  rcut_ = pRcut;
}

void TrackingAllPairs::UpdateTracking(bool pForceUpdate) {
  if (first_ || pForceUpdate) {
    for (int i = 0; i < nsimples_; ++i) {
      neighbors_[i].clear();
    }
    first_ = false;
    nupdates_++;
    // Loop over all particles, not double counting, and add
    // unique interactions to the mix
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

          neighbor_t new_neighbor;
          new_neighbor.idx_ = jdx;
          new_neighbor.rid_me_ = rid1;
          new_neighbor.rid_you_ = rid2;
          neighbors_[idx].push_back(new_neighbor);
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
}

void TrackingAllPairs::print() {
  // Is there really anything to print?
  printf("\tall pairs has nothing interesting to print\n");
}

void TrackingAllPairs::dump() {
  printf("\tall pairs has nothing interesting to dump\n");
}
