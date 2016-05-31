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
    n_updates_++;
    for (int i = 0; i < nparticles_; ++i) {
        neighbors_[i].clear();
    }
    
    // Update the cell list
    cell_list_.UpdateCellList(particles);
#ifdef DEBUG
    cell_list_.dump();
#endif
    
    // Call the cell update
    CellUpdate(particles);
#ifdef DEBUG
    dump();
#endif
    
    // Reset the total distance traveled for updates
    for (int i = 0; i < nparticles_; ++i) {
        memset((*particles)[i]->dr_tot, 0, sizeof(double)*3);
    }
}


// Attempt 2 of cell list update
// XXX: Work on stuff here when I want to
void
NeighborListCell::CellUpdate(std::vector<particle *> *particles) {

#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
        std::vector<int>* pid_to_cid = cell_list_.pidtocid();
        
        // Loop over all particles in this loop, and get the cell id
        // From that, only compare with particles that have a higher
        // pid than we do - this prevents double counting!
#if defined(_OPENMP)
#pragma omp for schedule(runtime) nowait
#endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Get our cell
            int cidx = (*pid_to_cid)[idx];
            auto cell1 = cell_list_[cidx];
            // Loop over other cells (including us) in the block of cells
            for (int cjdx = 0; cjdx < 27; ++cjdx) {
                auto cell2 = cell_list_[cell1->adj_cell_ids_[cjdx]];
                // Loop over it's particles
                for (int jdx = 0; jdx < cell2->nparticles_; ++jdx) {
                    int jjdx = cell2->idxlist_[jdx];
                    // ONLY DO THE CALCULATION IF THE OTHER PID IS HIGHER!!!
                    if (jjdx > idx) {
                        auto p1 = (*particles)[idx];
                        auto p2 = (*particles)[jjdx];
                        
                        double rx = buffmd::pbc(p1->x[0] - p2->x[0], boxby2_[0], box_[0]);
                        double ry = buffmd::pbc(p1->x[1] - p2->x[1], boxby2_[1], box_[1]);
                        double rz = buffmd::pbc(p1->x[2] - p2->x[2], boxby2_[2], box_[2]);
                        double rsq = rx*rx + ry*ry + rz*rz;
                        
                        if (rsq < rcs2_) {
                            neighbor_t new_neighbor;
                            new_neighbor.idx_ = jjdx;
                            neighbors_[idx].push_back(new_neighbor);
                        }
                    } // only do the calculation if the second id is higher
                } // cell2 particle list
            } // cells adjancent and equal to us
        } // pragma omp for reduction(+:epot) schedule(runtime) nowait
    } // pragma omp parallel
}


// Print neighbor list information
void
NeighborListCell::print() {
    // print out information like the cutoff radius, update distance
    // average number of neighbors
    printf("********\n");
    printf("Neighbor List Cells (threads:%d)\n", nthreads_);
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


// Dump everythign in gory detail
void
NeighborListCell::dump() {
    std::vector<int>* pid_to_cid = cell_list_.pidtocid();
    printf("\n********\n");
    printf("Dumping Neighbor List Cells!\n\n");
    printf("Neighbor List Cells -> {n:%d}{threads:%d}\n", nparticles_, nthreads_);
    // Loop through all particles
    printf("Printing neighbors\n");
    for (int idx = 0; idx < nparticles_; ++idx) {
        int cidx = (*pid_to_cid)[idx];
        printf("\tp(%d) -> {n:%lu}{cell:%d}\n", idx, neighbors_[idx].size(), cidx);
        printf("\t neighbors: ");
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
            int jdx = nldx->idx_;
            printf("%d,", jdx);
        }
        printf("\n");
    }
    printf("NOW DUMPING SUB CELL LIST!\n");
    cell_list_.dump();
    printf("DONE DUMPING SUB CELL LIST!\n");
    printf("\n********\n");
}
