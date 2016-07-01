// implementation for neighbor lists (all pairs)

#include <chrono>
#include <climits>
#include <set>

#include "neighbor_list_ap.h"

#include "minimum_distance.h"

void
NeighborListAP::CreateSubstructure(double pRcut) {
  auto start = std::chrono::steady_clock::now();

  rcut_ = pRcut;
  rcut2_ = rcut_*rcut_;

  // Make sure that if we don't have a skin to make it
  // something...
  if (skin_ <= 0.0)
    skin_ = 1.0;

  skin2_ = skin_*skin_;
  rcs2_ = rcut2_ + skin2_ + 2*rcut_*skin_;
  half_skin2_ = 0.25*skin2_;

  neighbors_ = new nl_list[nparticles_];

  for (int idim = 0; idim < ndim_; ++idim) {
    boxby2_[idim] = 0.5*box_[idim];
  }

  n_updates_ = 0;

  auto end = std::chrono::steady_clock::now();
  std::cout << "NeighborListAP::CreateSubstructure: " << std::chrono::duration<double, std::milli> (end-start).count() << "ms\n";
}


// Check to see if the neighbor list needs updating based on
// accumulators
void
NeighborListAP::CheckNeighborList(bool pForceUpdate) {
  nl_update_ = pForceUpdate;
  for (auto spec = species_->begin(); spec!=species_->end(); ++spec) {
    if (nl_update_)
    break;
    double dr2 = (*spec)->GetDrMax();
    nl_update_ = dr2 > half_skin2_;
  }
  if (nl_update_) {
    UpdateNeighborList();
  }
}


// Update the neighbor list
void
NeighborListAP::UpdateNeighborList() {
  // Purge the neighbor list
  n_updates_++;
  for (int i = 0; i < nparticles_; ++i) {
    neighbors_[i].clear();
  }

  // Call the all pairs update
  AllPairsUpdate();

  // Reset the accumulators
  for (auto spec = species_->begin(); spec!=species_->end(); ++spec) {
    (*spec)->ZeroDr();
  }
}


// All pairs update routine
void
NeighborListAP::AllPairsUpdate() {
  // Loop over all pairs to build the neighbor list (simple, but O(N^2))
  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    std::set< std::pair<int, int> > rid_interactions;
    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int idx = 0; idx < nparticles_ - 1; ++idx) {
      for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
        auto p1 = (*simples_)[idx];
        auto p2 = (*simples_)[jdx];
        int rid1 = p1->GetRID();
        int rid2 = p2->GetRID();
        // Only insert object if rigid ID is unique
        // std::set only allows insertion of unique elements
        // Check if insertion fails because RID is not unique
        if (!rid_interactions.insert(std::make_pair(rid1,rid2)).second || 
            !rid_interactions.insert(std::make_pair(rid2,rid1)).second)
          continue;
        if (p1->GetCID() == p2->GetCID()) continue;

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
    } // pragma omp for schedule(runtime) nowait
  } // pragma omp parallel
  #ifdef ENABLE_OPENMP
  {
    // Clean up any existing rid-rid interactions
    std::set< std::pair<int, int> > rid_interactions;
    for (int idx=0; idx<nparticles_; ++idx) {
      for (auto nb = neighbors_[idx].begin(); nb != neighbors_[idx].end();) {
        if (!rid_interactions.insert(std::make_pair(nb->rid_me_,nb->rid_you_)).second ||
            !rid_interactions.insert(std::make_pair(nb->rid_you_,nb->rid_me_)).second) {
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


// print
void
NeighborListAP::print() {
  printf("********\n");
  printf("%s ->\n", name_.c_str());
  printf("\t{rcut: %2.2f}, {skin: %2.2f} = {rcs2: %2.2f}, {half_skin2:%2.2f}\n", rcut_, skin_, rcs2_, half_skin2_);
  int ntotlist = 0, maxlist = 0, minlist = INT_MAX;
  for (int i = 0; i < nparticles_ -1; ++i) {
    ntotlist += neighbors_[i].size();
    maxlist = std::max(maxlist, (int)neighbors_[i].size());
    minlist = std::min(minlist, (int)neighbors_[i].size());
  }
  printf("\tStats: {min: %d}, {max: %d}, {avg: %2.2f}\n", minlist, maxlist, (float)ntotlist/(float)nparticles_);
}


// dump gory details
void
NeighborListAP::dump() {
  #ifdef DEBUG
  printf("********\n");
  printf("%s -> dump\n", name_.c_str());
  for (int idx = 0; idx < nparticles_; ++idx) {
    printf("\t[%d] -> [", idx);
    for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
      int jdx = nldx->idx_;
      printf("%d,", jdx);
    }
    printf("]\n");
  }
  #endif
}
