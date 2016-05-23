// Implementation for the neighbor list base ideas

#include <cstdlib>
#include <climits>
#include <algorithm>

#include "neighbor_list_cell.h"

NeighborListCell::~NeighborListCell() {
    //delete(neighbors_);
}


// Create the neighbor list
void
NeighborListCell::CreateNeighborList(int pN, double pRcut, double pSkin, double pBox[3]) {
    // Call the base constructor
    NeighborList::CreateNeighborList(pN, pRcut, pSkin, pBox);
    
    // Create the cell list
    cell_list_.CreateCellList(pN, pRcut, pSkin, pBox);
}


// Update the neighbor list
void
NeighborListCell::UpdateNeighborList(std::vector<particle *> *particles) {
    // Purge the neighbor list
    printf("Starting UpdateNeighborList\n");
    n_updates_++;
    for (int i = 0; i < nparticles_; ++i) {
        neighbors_[i].clear();
    }
    
    // Update the cell list
    cell_list_.UpdateCellList(particles);
    
    // Call the cell update
    CellUpdate(particles);
    
    // Reset the total distance traveled for updates
    for (int i = 0; i < nparticles_; ++i) {
        memset((*particles)[i]->dr_tot, 0, sizeof(double)*3);
    }
    printf("Finished UpdateNeighborList\n");
}


// Cell update update
// Doesn't work, some sort of synchronizationissue with the neighbor list
/*void
NeighborListCell::CellUpdate(std::vector<particle *> *particles) {
    // Loop over the cells to build the neighbor list
    printf("Starting CellUpdate\n");
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
        int tid;
#if defined(_OPENMP)
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        
        // Check within my own cell
        int ncells = cell_list_.ncells();
        printf("Checking within my own cell, tid:%d\n", tid);
        for (int cidx = 0; cidx < ncells; cidx += nthreads_) {
            // set the index
            int cjdx = cidx + tid;
            if (cjdx >= ncells) break;
            
            // Get the actual cell
            auto c1 = cell_list_[cjdx];
            // Loop over particles in said cell
            for (int pidx1 = 0; pidx1 < c1->nparticles_ - 1; ++pidx1) {
                int idx = c1->idxlist_[pidx1];
                auto p1 = (*particles)[idx];
                
                // Get my interacting partner
                for(int pidx2 = pidx1 + 1; pidx2 < c1->nparticles_; ++pidx2) {
                    int jdx = c1->idxlist_[pidx2];
                    auto p2 = (*particles)[jdx];
                    
                    //printf("{c1:%d,c2:SAME}->{p1:%d,p2:%d}\n", c1->cell_id_, idx, jdx);
                    
                    double rx = buffmd::pbc(p1->x[0] - p2->x[0], boxby2_[0], box_[0]);
                    double ry = buffmd::pbc(p1->x[1] - p2->x[1], boxby2_[1], box_[1]);
                    double rz = buffmd::pbc(p1->x[2] - p2->x[2], boxby2_[2], box_[2]);
                    double rsq = rx*rx + ry*ry + rz*rz;
                    
                    if (rsq < rcs2_) {
                        //printf("\t->Adding to neighbor list n[%d] = %d\n", idx, jdx);
                        neighbor_t new_neighbor;
                        new_neighbor.idx_ = jdx;
                        neighbors_[idx].push_back(new_neighbor);
                    }
                } // check interaction partner particles
            } // check the actual particles
        }  // Cell loop
        
        // Interactions across different cells
        printf("Checking with interacting cells, tid:%d\n", tid);
        int npairs = cell_list_.npairs();
        for (int pairidx = 0; pairidx < npairs; pairidx += nthreads_) {
            int pairjdx = pairidx + tid;
            if (pairjdx >= npairs) break;
            auto cell1 = cell_list_[cell_list_.plist(2*pairjdx  )];
            auto cell2 = cell_list_[cell_list_.plist(2*pairjdx+1)];
            
            for (int pidx1 = 0; pidx1 < cell1->nparticles_; ++pidx1) {
                int idx = cell1->idxlist_[pidx1];
                auto p1 = (*particles)[idx];
                
                for (int pidx2 = 0; pidx2 < cell2->nparticles_; ++pidx2) {
                    int jdx = cell2->idxlist_[pidx2];
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
                } // Second particle
            } // First particle
        } // Pairs loops
#if defined(_OPENMP)
#pragma omp barrier
#endif
    } // omp loop
    printf("Finished CellUpdate\n");
}*/


// Attempt 2 of cell list update
// XXX: Work on stuff here when I want to
void
NeighborListCell::CellUpdate(std::vector<particle *> *particles) {
    printf("Starting CellUpdate\n");
    printf("WARNING: NOT IMPLEMENTED RIGHT NOW, EXITING!\n");
    exit(1);
    
    for (int idx = 0; idx < nparticles_; ++idx) {
        // Get the cell containing this particle
        //int cidx = (*particles)[idx]->cellid;
    }
    
    printf("Finished CellUpdate\n");
}


// Print neighbor list information
void
NeighborListCell::print() {
    // print out information like the cutoff radius, update distance
    // average number of neighbors
    printf("********\n");
    printf("Neighbor List Cells\n");
    printf("\t{rcut:%2.2f}, {skin:%2.2f}, {rcs2:%2.2f}\n", rcut_, skin_, rcs2_);
    int ntotlist = 0, maxlist = 0, minlist = INT_MAX;
    // The -1 is important because the last neighbor list is empty when we aren't
    // keeping every i -> j -> i pair (only keeping half)
    for (int i = 0; i < nparticles_ - 1; ++i) {
        ntotlist += neighbors_[i].size();
        maxlist = std::max(maxlist, (int)neighbors_[i].size());
        minlist = std::min(minlist, (int)neighbors_[i].size());
    }
    printf("\t{Neighbors stats -> {min:%d}, {max:%d}, {avg:%2.2f}\n", minlist, maxlist, (float)ntotlist/(float)nparticles_);
    printf("\t{NeighborList size: %.1fMb}\n", (float)GetMemoryFootprint()/1024/1024.);
    // Print out the interesting information from the cell list
    cell_list_.CheckCellList();
}
