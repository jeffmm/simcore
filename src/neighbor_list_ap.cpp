// implementation for neighbor lists (all pairs)

#include <chrono>

#include "neighbor_list_ap.h"

#include "minimum_distance.h"

void
NeighborListAP::CreateSubstructure(double pRcut) {
    auto start = std::chrono::steady_clock::now();

    rcut_ = pRcut;
    rcut2_ = rcut_*rcut_;

    // XXX: CJE for now overwrite skin to something until we add
    // it via params
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
    printf("-> {rcut: %2.2f}, {skin: %2.2f} = {rcs: %2.2f}, {half_skin:%2.2f}\n", rcut_, skin_, rcs2_, half_skin2_);

    std::cout << "NeighborListAP::CreateSubstructure: " << std::chrono::duration<double, std::milli> (end-start).count() << "ms\n";
}


// Check to see if the neighbor list needs updating based on
// accumulators
void
NeighborListAP::CheckNeighborList(bool pForceUpdate) {
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
NeighborListAP::UpdateNeighborList() {
    // Purge the neighbor list
    n_updates_++;
    for (int i = 0; i < nparticles_; ++i) {
        neighbors_[i].clear();
    }

    // Call the all pairs update
    AllPairsUpdate();

    // Reset the accumulators
    for (int i = 0; i < nparticles_; ++i) {
        auto part = simples_[i];
        part->ZeroDrTot();
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
        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
                auto p1 = simples_[idx];
                auto p2 = simples_[jdx];

                // Minimum distance
                interactionmindist idm;
                MinimumDistance(p1, p2, idm, ndim_, nperiodic_, space_);

                if (idm.dr_mag2 < rcs2_) {
                    neighbor_t new_neighbor;
                    new_neighbor.idx_ = jdx;
                    neighbors_[idx].push_back(new_neighbor);
                }
            }
        } // pragma omp for schedule(runtime) nowait
    }
}
