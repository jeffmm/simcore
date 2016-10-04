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
    //std::cout << "Creating unique rids\n";
    unique_rids_->clear();
    for (int i = 0; i < nsimples_; ++i) {
      neighbors_[i].clear();
      // RID checker
      unique_rids_->insert((*simples_)[i]->GetRID());
    }
    nrigids_ = (int)unique_rids_->size();
    //std::cout << "Found " << nrigids_ << " rigids\n";
    maxrigid_ = *(unique_rids_->rbegin());
    //std::cout << "Max rigid rid: " << maxrigid_ << std::endl << std::flush;

    first_ = false;
    nupdates_++;

    // Loop over particles, 2 way nl, and add unique inter
    // actions
    rid_self_check_->clear();
    //rid_self_check_->resize(nrigids_);
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

        //std::cout << "Checking simple: " << p1->GetOID() << ", idx: " << idx << ", rid: " << rid1 << std::endl;
        //p1->Dump();

        //std::cout << "rid self check: ";
        //for (auto ridit = rid_self_check_->begin(); ridit != rid_self_check_->end(); ++ridit) {
        //  std::cout << *ridit << ", ";
        //}
        //std::cout << "done" << std::endl << std::flush;

        bool should_exit = false;
        #ifdef ENABLE_OPENMP
        #pragma omp critical
        #endif
        {
          should_exit = (*rid_self_check_)[rid1];
          //std::cout << "Checking rid: " << rid1 << " gives " << (should_exit ? "true" : "false") << std::endl << std::flush;
          (*rid_self_check_)[rid1] = true;
        }
        if (should_exit) {
          //std::cout << "  Already checked rid " << rid1 << ", continuing\n";
          continue;
        }

        // Get the threadlocal rid_check_local_
        //std::cout << "--Getting threadlocal tid " << tid << std::endl << std::flush;
        std::unordered_set<int>* rid_check_local = rid_check_local_[tid];
        //std::cout << "--Clearing threadlocal: " << rid_check_local << std::endl << std::flush;
        rid_check_local->clear();
        //std::cout << "--Cleared threadlocal\n" << std::flush;

        for (int jdx = 0; jdx < nsimples_; ++jdx) {
          if (idx == jdx) continue;
          auto p2 = (*simples_)[jdx];
          int rid2 = p2->GetRID();
          if (rid1 == rid2) continue;

          // Check if there is even an interaction
          PotentialBase *pot1 = potentials_->GetPotentialExternal(p1->GetSID(), p2->GetSID());
          PotentialBase *pot2 = potentials_->GetPotentialInternal(p1->GetOID(), p2->GetOID());
          if ((pot1 == nullptr) && (pot2 == nullptr)) continue;

          //std::cout << "    Checking simple: " << p2->GetOID() << ", jdx: " << jdx << ", rid: " << rid2 << std::endl << std::flush;
          if (rid_check_local->count(rid2)) {
            //std::cout << "    Did find rid " << rid2 << ", continuing\n" << std::flush;
            continue;
          }
          //std::cout << "    Did not find rid " << rid2 << std::endl << std::flush;
          rid_check_local->insert(rid2);

          //std::cout << "    blah: " << std::flush;
          for (auto blah = rid_check_local->begin(); blah != rid_check_local->end(); ++blah) {
            //std::cout << *blah << ", " << std::flush;
          }
          //std::cout << "blahend\n" << std::flush;

          // Don't do a distance check, just add
          neighbor_t new_neighbor;
          new_neighbor.idx_ = jdx;
          new_neighbor.rid_me_ = rid1;
          new_neighbor.rid_you_ = rid2;
          //std::cout << "    Creating new neighbor idx: " << new_neighbor.idx_ << ", rid1: " << new_neighbor.rid_me_
          //  << ", rid2: " << new_neighbor.rid_you_ << std::endl << std::flush;
          neighbors_[idx].push_back(new_neighbor);
        }
      } // pragma omp for
    }
  }

  //std::cout << "Finished update tracking\n" <<std::flush;
}

/*
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
}
*/

void TrackingAllPairs::print() {
  // Is there really anything to print?
  printf("\tall pairs has nothing interesting to print\n");
}

void TrackingAllPairs::dump() {
  printf("\tall pairs has nothing interesting to dump\n");
}
