// Implementation for the neighbor list base ideas

#include <cstdlib>
#include <climits>
#include <algorithm>

#include "neighbor_list.h"

NeighborList::~NeighborList() {
    delete(neighbors_);
}


// Create the neighbor list
void
NeighborList::CreateNeighborList(int pN, double pRcut, double pSkin, double pBox[3]) {
    nparticles_ = pN;
    rcut_ = pRcut;
    rcut2_ = rcut_*rcut_;
    skin_ = pSkin;
    skin2_ = skin_*skin_;
    rcs2_ = rcut2_ + skin2_ + 2*rcut_+skin_;
    half_skin2_ = 0.25*skin2_;
    memcpy(box_, pBox, 3*sizeof(double));
    
    neighbors_ = new nl_list[nparticles_];
    
    for (int i = 0; i < 3; ++i) {
        boxby2_[i] = 0.5 * box_[i];
    }
    
    n_updates_ = 0;
}


// Check the conditions if the neighbor list needs updating!
void
NeighborList::CheckNeighborList(std::vector<particle*>* particles) {
    for (int idx = 0; idx < nparticles_; ++idx) {
        auto p1 = (*particles)[idx];
        double dr2 = p1->dr_tot[0]*p1->dr_tot[0] +
                     p1->dr_tot[1]*p1->dr_tot[1] +
                     p1->dr_tot[2]*p1->dr_tot[2];
        nl_update_ = dr2 > half_skin2_;
        if (nl_update_) {
            break;
        }
    }
    
    if (nl_update_) {
        UpdateNeighborList(particles);
    }
}


// Update the neighbor list
void
NeighborList::UpdateNeighborList(std::vector<particle *> *particles) {
    // Purge the neighbor list
    n_updates_++;
    for (int i = 0; i < nparticles_; ++i) {
        neighbors_[i].clear();
    }
    
    // Call the all pairs update
    AllPairsUpdate(particles);
    
    // Reset the total distance traveled for updates
    for (int i = 0; i < nparticles_; ++i) {
        memset((*particles)[i]->dr_tot, 0, sizeof(double)*3);
    }
}


// All pairs update
void
NeighborList::AllPairsUpdate(std::vector<particle *> *particles) {
    // Loop over all pairs and build the neighbor list
#pragma omp parallel
    {
#pragma omp for schedule(runtime) nowait
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
                auto p1 = (*particles)[idx];
                auto p2 = (*particles)[jdx];
                
                double rx = buffmd::pbc(p1->x[0] - p2->x[0], boxby2_[0], box_[0]);
                double ry = buffmd::pbc(p1->x[1] - p2->x[1], boxby2_[1], box_[1]);
                double rz = buffmd::pbc(p1->x[2] - p2->x[2], boxby2_[2], box_[2]);
                double rsq = rx*rx + ry*ry + rz*rz;
                
                if (rsq < rcs2_) {
                    neighbor_t new_neighbor;
                    new_neighbor.idx_ = jdx;
                    neighbors_[idx].push_back(new_neighbor);
                }
            }
        }
    }
}


// Print neighbor list information
void
NeighborList::print() {
    // print out information like the cutoff radius, update distance
    // average number of neighbors
    printf("********\n");
    printf("Neighbor List All Pairs\n");
    printf("\t{rcut:%.2f}, {skin:%.2f}, {rcs2:%.2f}\n", rcut_, skin_, rcs2_);
    int ntotlist = 0, maxlist = 0, minlist = INT_MAX;
    // The -1 is important because the last neighbor list is empty when we aren't
    // keeping every i -> j -> i pair (only keeping half)
    for (int i = 0; i < nparticles_ - 1; ++i) {
        ntotlist += neighbors_[i].size();
        maxlist = std::max(maxlist, (int)neighbors_[i].size());
        minlist = std::min(minlist, (int)neighbors_[i].size());
    }
    printf("\t{Neighbors stats -> {min:%d}, {max:%d}, {avg:%.2f}\n", minlist, maxlist, (float)ntotlist/(float)nparticles_);
    printf("\t{NeighborList size: %.1fMb}\n", (float)GetMemoryFootprint()/1024/1024.);
}


// Memory footprint stuffs
unsigned long
NeighborList::GetMemoryFootprint() {
    auto mysize = sizeof(*this);
    unsigned long subsize = 0;
    for (int i = 0; i < nparticles_; ++i) {
        subsize += neighbors_[i].size();
    }
    subsize *= sizeof(neighbor_t);
    return mysize + subsize;
}
