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
    std::cout << "Creating unique rids\n";
    std::unordered_set<int> unique_rids;
    for (int i = 0; i < nsimples_; ++i) {
      neighbors_[i].clear();
      // RID checker
      unique_rids.insert((*simples_)[i]->GetRID());
    }
    nrigids_ = (int)unique_rids.size();
    std::cout << "Found " << nrigids_ << " rigids\n";

    first_ = false;
    nupdates_++;

    // Loop over particles, 2 way nl, and add unique inter
    // actions
    std::vector<bool> rid_self_check(nrigids_, false);
    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef ENABLE_OPENMP
      #pragma omp for schedule(runtime) nowait
      #endif
      for (int idx = 0; idx < nsimples_; ++idx) {
        std::vector<int> rid_interactions_local;
        auto p1 = (*simples_)[idx];
        int rid1 = p1->GetRID();

        std::cout << "Checking simple: " << p1->GetOID() << ", idx: " << idx << ", rid: " << rid1 << std::endl;

        bool should_exit = false;
        #ifdef ENABLE_OPENMP
        #pragma omp critical
        #endif
        {
          should_exit = rid_self_check[rid1];
          rid_self_check[rid1] = true;
        }
        if (should_exit) {
          std::cout << "  Already checked rid " << rid1 << ", continuing\n";
          continue;
        }

        for (int jdx = 0; jdx < nsimples_; ++jdx) {
          if (idx == jdx) continue;
          auto p2 = (*simples_)[jdx];
          int rid2 = p2->GetRID();
          if (rid1 == rid2) continue;
          std::cout << "    Checking simple: " << p2->GetOID() << ", jdx: " << jdx << ", rid: " << rid2 << std::endl;
          auto ridit = find(rid_interactions_local.begin(), rid_interactions_local.end(), rid2);
          if (ridit != rid_interactions_local.end()) {
            std::cout << "    Found that rid " << *ridit << " already exists, continuing\n";
            continue;
          }
          std::cout << "    Did not find rid " << rid2 << std::endl << std::flush;
          rid_interactions_local.push_back(rid2);

          std::cout << "blah: " << std::flush;
          for (auto blah = rid_interactions_local.begin(); blah != rid_interactions_local.end(); ++blah) {
            std::cout << *blah << ", " << std::flush;
          }
          std::cout << "blahend\n" << std::flush;

          // Don't do a distance check, just add
          neighbor_t new_neighbor;
          new_neighbor.idx_ = jdx;
          new_neighbor.rid_me_ = rid1;
          new_neighbor.rid_you_ = rid2;
          std::cout << "    Creating new neighbor idx: " << new_neighbor.idx_ << ", rid1: " << new_neighbor.rid_me_
            << ", rid2: " << new_neighbor.rid_you_ << std::endl << std::flush;
          neighbors_[idx].push_back(new_neighbor);
        }
      } // pragma omp for
    }
  }
}

/*void TrackingAllPairs::UpdateTracking(bool pForceUpdate) {
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
}*/

void TrackingAllPairs::print() {
  // Is there really anything to print?
  printf("\tall pairs has nothing interesting to print\n");
}

void TrackingAllPairs::dump() {
  printf("\tall pairs has nothing interesting to dump\n");
}
