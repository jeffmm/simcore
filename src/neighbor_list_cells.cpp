// implementation for neighbor lists (all pairs)

#include <chrono>
#include <climits>

#include "neighbor_list_cells.h"

#include "minimum_distance.h"

// Init must call the cell list stuff too
void
NeighborListCells::Init(space_struct *pSpace, double pSkin) {
    // Call the base init
    printf("NeighborListCells::Init begin\n");
    ForceSubstructureBase::Init(pSpace, pSkin);

    cell_list_.Init(pSpace, pSkin);
    printf("NeighborListCells::Init end\n");
}


// Load flat simples must be overridden
void
NeighborListCells::LoadFlatSimples(std::vector<Simple*> pSimples) {
    printf("NeighborListCells::LoadFlatSimples begin\n");
    // Call the base class version
    ForceSubstructureBase::LoadFlatSimples(pSimples);

    // Call the cell list version
    cell_list_.LoadFlatSimples(pSimples);
    printf("NeighborListCells::LoadFlatSimples end\n");
}

void
NeighborListCells::CreateSubstructure(double pRcut) {
    printf("NeighborListCells::CreateSubstructure begin\n");
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

    // Also create the cell list substructure here!!!!!!!
    cell_list_.CreateSubstructure(pRcut);
    cell_list_.UpdateCellList();

    std::cout << "NeighborListCells::CreateSubstructure: " << std::chrono::duration<double, std::milli> (end-start).count() << "ms\n";
    printf("NeighborListCells::CreateSubstructure end\n");
}


// Check to see if the neighbor list needs updating based on
// accumulators
void
NeighborListCells::CheckNeighborList(bool pForceUpdate) {
    nl_update_ = pForceUpdate;
    for (int idx = 0; idx < nparticles_; ++idx) {
        if (nl_update_)
            break;

        auto part = simples_[idx];
        auto pr = part->GetDrTot();

        double dr2 = 0.0;
        for (int i = 0; i < ndim_; ++i) {
            dr2 += pr[i]*pr[i];
        }
        nl_update_ = dr2 > half_skin2_;
    }

    if (nl_update_) {
        UpdateNeighborList();
    }
}


// Update the neighbor list
void
NeighborListCells::UpdateNeighborList() {
    // Purge the neighbor list
    n_updates_++;
    for (int i = 0; i < nparticles_; ++i) {
        neighbors_[i].clear();
    }

    // Update the cell list
    cell_list_.UpdateCellList();

    // Call the all pairs update
    CellsUpdate();

    // Reset the accumulators
    for (int i = 0; i < nparticles_; ++i) {
        auto part = simples_[i];
        part->ZeroDrTot();
    }
}


// Cell list based update routine
void
NeighborListCells::CellsUpdate() {
    // Loop over cells to build neighbor list (O(N))
    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<int> *pid_to_cid = cell_list_.pidtocid();

        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            // Get our cell
            int cidx = (*pid_to_cid)[idx];
            auto cell1 = cell_list_[cidx];
            // Loop over all other cells (including this one)
            int nadj = cell_list_.nadj();
            for (int cjdx = 0; cjdx < nadj; ++cjdx) {
                auto cell2 = cell_list_[cjdx];
                // Loop over cell 2 particles
                for (int jdx = 0; jdx < cell2->nparticles_; ++jdx) {
                    int jjdx = cell2->idxlist_[jdx];
                    if (jjdx > idx) {
                        auto part1 = simples_[idx];
                        auto part2 = simples_[jjdx];

                        // Minimum distance
                        interactionmindist idm;
                        MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);

                        if (idm.dr_mag2 < rcs2_) {
                            neighbor_t new_neighbor;
                            new_neighbor.idx_ = jjdx;
                            neighbors_[idx].push_back(new_neighbor);
                        }
                    }
                }
            }
        } // pragma omp for schedule(runtime) nowait
    } // pragma omp parallel
}


// print
void
NeighborListCells::print() {
    printf("********\n");
    printf("%s ->\n", name_.c_str());
    printf("\t{rcut: %2.2f}, {skin: %2.2f} = {rcs2: %2.2f}, {half_skin2:%2.2f}\n", rcut_, skin_, rcs2_, half_skin2_);
    int ntotlist = 0, maxlist = 0, minlist = INT_MAX;
    for (int i = 0; i < nparticles_ -1; ++i) {
        ntotlist += neighbors_[i].size();
        maxlist = std::max(maxlist, (int)neighbors_[i].size());
        minlist = std::min(minlist, (int)neighbors_[i].size());
    }
    printf("\t--------\n");
    cell_list_.print();
    printf("\t--------\n");
    printf("\tStats: {min: %d}, {max: %d}, {avg: %2.2f}\n", minlist, maxlist, (float)ntotlist/(float)nparticles_);
}


// dump gory details
void
NeighborListCells::dump() {
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
    cell_list_.dump();
    #endif
}
