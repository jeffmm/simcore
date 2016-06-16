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
                MinimumDistance(p1, p2, idm);

                if (idm.dr_mag2 < rcs2_) {
                    neighbor_t new_neighbor;
                    new_neighbor.idx_ = jdx;
                    neighbors_[idx].push_back(new_neighbor);
                }
            }
        } // pragma omp for schedule(runtime) nowait
    }
}

// Find the minimum distance beween two particles
void NeighborListAP::MinimumDistance(Simple* o1, Simple* o2, interactionmindist& imd) {
    double const * const r1 = o1->GetPosition();
    double const * const s1 = o1->GetScaledPosition();
    double const * const u1 = o1->GetOrientation();
    double const * const r2 = o2->GetPosition();
    double const * const s2 = o2->GetScaledPosition();
    double const * const u2 = o2->GetOrientation();
    double const l1 = o1->GetLength();
    double const l2 = o2->GetLength();
    double const d1 = o1->GetDiameter();
    double const d2 = o2->GetDiameter();
    /* TODO: Think about how best to do this for general shapes, like 2d
       polygons that can represent the local surface of more complex 3d
       shapes. Perhaps assume all local surface to be triangular polygons.*/
    imd.dr_mag2 = 0;
    std::fill(imd.dr, imd.dr+3, 0.0);
    std::fill(imd.contact1, imd.contact1+3, 0.0);
    std::fill(imd.contact2, imd.contact2+3, 0.0);
    imd.buffer_mag = 0.5*(d1+d2);
    imd.buffer_mag2 = imd.buffer_mag*imd.buffer_mag;
    if (l1 == 0 && l2 == 0)
        min_distance_point_point(ndim_, nperiodic_, space_->unit_cell, 
                                 r1, s1, r2, s2, imd.dr, &imd.dr_mag2);
    else if (l1 == 0 && l2 > 0) 
        min_distance_sphere_sphero(ndim_, nperiodic_, space_->unit_cell,
                                   r1, s1, r2, s2, u2, l2,
                                   imd.dr, &imd.dr_mag2, imd.contact2);
    else if (l1 > 0 && l2 == 0) 
        min_distance_sphere_sphero(ndim_, nperiodic_, space_->unit_cell,
                                   r2, s2, r1, s1, u1, l1,
                                   imd.dr, &imd.dr_mag2, imd.contact1);
    else if (l1 > 0 && l2 > 0)
        min_distance_sphero(ndim_, nperiodic_, space_->unit_cell,
                            r1, s1, u1, l1, r2, s2, u2, l2,
                            imd.dr, &imd.dr_mag2, imd.contact1, imd.contact2);
    imd.dr_mag = sqrt(imd.dr_mag2);
}
