// Implementation for force microcells

#include <cassert>

#include "force_neighborlistap.h"

// Overridden init method to call initmp
void
ForceNeighborListAP::Init(space_struct* pSpace, double pSkin) {
    // Override this to call base class, then initmp
    ForceBase::Init(pSpace, pSkin);
    InitMP();
}


// Overridden load simples method (needed for finalize)
void
ForceNeighborListAP::LoadSimples(std::vector<SpeciesBase*> pSpecies) {
    // Run the base class version
    ForceBase::LoadSimples(pSpecies);

    // Load the flat simples into the cell list
    neighbor_list_.LoadFlatSimples(simples_);
}

void
ForceNeighborListAP::InitMP() {
    neighbor_list_.Init(space_, skin_);
}


void
ForceNeighborListAP::Finalize() {
    // Should be ready to finalize if we have loaded and init'ed
    neighbor_list_.CreateSubstructure(max_rcut_); 
    //UpdateScheme();
    neighbor_list_.CheckNeighborList(true);
    initialized_ = true;
}


void
ForceNeighborListAP::UpdateScheme() {
    neighbor_list_.CheckNeighborList();
}


void
ForceNeighborListAP::Interact() {
    if (!initialized_) {
        printf("ERROR: Finalized was not run for ForceNeighborListAPs, exiting!\n");
        exit(1);
    }


    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
        int tid;
        double **fr = new double*[3];
        double **tr = new double*[3];
        double *pr_energy;

        auto neighbors = neighbor_list_.GetNeighbors();

        #ifdef ENABLE_OPENMP
        tid = omp_get_thread_num();
        #else
        tid = 0;
        #endif

        // Set up the pointers to the force  and torque superarrays
        // Do all 3 dimensions just to be safe
        for (int i = 0; i < 3; ++i) {
            fr[i] = frc_ + ((3*tid+i)*nparticles_);
            tr[i] = trqc_ + ((3*tid+i)*nparticles_);
        }

        for (int i = 0; i < 3*nparticles_; ++i) {
            (*fr)[i] = 0.0;
            (*tr)[i] = 0.0;
        }

        pr_energy = prc_energy_ + (tid * nparticles_);
        for (int i = 0; i < nparticles_; ++i) {
            pr_energy[i] = 0.0;
        }

        assert(nparticles_ == simples_.size());
        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Iterate over our neighbors
            for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
                int jdx = nldx->idx_;
                auto part1 = simples_[idx];
                auto part2 = simples_[jdx];

                // Exclude composite object interactions
                if (part1->GetCID() == part2->GetCID()) continue;

                // Calculate the potential here
                PotentialBase *pot = potentials_.GetPotential(part1->GetSID(), part2->GetSID());
                if (pot == nullptr) continue;
                // Minimum distance here@@@@!!!!
                // XXX: CJE ewwwwwww, more elegant way?
                interactionmindist idm;
                MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);
                if (idm.dr_mag2 > pot->GetRCut2()) continue;

                // Obtain the mapping between particle oid and position in the force superarray
                auto oid1x = oid_position_map_[part1->GetOID()];
                auto oid2x = oid_position_map_[part2->GetOID()];

                // Fire off the potential calculation
                double fepot[4];
                pot->CalcPotential(idm.dr, idm.dr_mag, idm.buffer_mag, fepot);

                // Do the potential energies
                pr_energy[oid1x] += fepot[ndim_];
                pr_energy[oid2x] += fepot[ndim_];

                // Do the forces
                for (int i = 0; i < ndim_; ++i) {
                    fr[i][oid1x] += fepot[i];
                    fr[i][oid2x] -= fepot[i];
                }

                // Calculate the torques
                double tau[3];
                cross_product(idm.contact1, fepot, tau, ndim_);
                for (int i = 0; i < ndim_; ++i) {
                    tr[i][oid1x] += tau[i];
                }
                cross_product(idm.contact2, fepot, tau, ndim_);
                for (int i = 0; i < ndim_; ++i) {
                    tr[i][oid2x] -= tau[i];
                }
            }
        } // pragma omp for schedule(runtime) nowait

        // Reduce once all threads have finished
        #ifdef ENABLE_OPENMP
        #pragma omp barrier
        #endif
        // Everything is doubled to do energy
        // XXX: CJE maybe fix this later?
        int i = 1 + (3 * nparticles_ / nthreads_);
        int ii = 1 + (nparticles_ / nthreads_);

        int fromidx = tid * i;
        int fromidx2 = tid * ii;

        int toidx = fromidx + i;
        int toidx2 = fromidx2 + ii;

        if (toidx > 3*nparticles_) toidx = 3*nparticles_;
        if (toidx2 > nparticles_) toidx2 = nparticles_;

        // Reduce the forces
        for (i = 1; i < nthreads_; ++i) {
            int offs;
            offs = 3 * i * nparticles_;

            for (int j = fromidx; j < toidx; ++j) {
                frc_[j] += frc_[offs+j];
                trqc_[j] += trqc_[offs+j];
            }
        }

        // Reduce the energies
        for (ii = 1; ii < nthreads_; ++ii) {
            int offs;
            offs = ii * nparticles_;

            for (int jj = fromidx2; jj < toidx2; ++jj) {
                prc_energy_[jj] += prc_energy_[offs + jj];
            }
        }

        //free(fr);
        //free(tr);
        delete[] fr;
        delete[] tr;
    } // pragma omp parallel

    // XXX: CJE this is where we tell the particles what to do?
    for (int i = 0; i < nparticles_; ++i) {
        auto part = simples_[i];
        int oidx = oid_position_map_[part->GetOID()];
        double subforce[3] = {0.0, 0.0, 0.0};
        double subtorque[3] = {0.0, 0.0, 0.0};
        for (int idim = 0; idim < ndim_; ++idim) {
            subforce[idim] = frc_[idim*nparticles_+oidx];
            subtorque[idim] = trqc_[idim*nparticles_+oidx]; 
        }
        part->AddForceTorqueEnergy(subforce, subtorque, prc_energy_[oidx]);
    }
}
