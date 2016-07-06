// Implementation for force neighbor list all pairs

#include <cassert>

#include "force_neighborlist_cells.h"

// Overridden init method to call initmp
void
ForceNeighborListCells::Init(space_struct* pSpace, std::vector<SpeciesBase*> *pSpecies, std::vector<Simple*> *pSimples, double pSkin) {
    // Override this to call base class, then initmp
    ForceBase::Init(pSpace, pSpecies, pSimples, pSkin);
    InitMP();
}


// Overridden load simples method (needed for finalize)
void
ForceNeighborListCells::LoadSimples() {
    // Run the base class version
    ForceBase::LoadSimples();

    // Load the flat simples into the cell list
    neighbor_list_.LoadFlatSimples();
}

void
ForceNeighborListCells::InitMP() {
    neighbor_list_.Init(space_, species_, simples_, skin_);
}


void
ForceNeighborListCells::Finalize() {
    // Should be ready to finalize if we have loaded and init'ed
    neighbor_list_.CreateSubstructure(max_rcut_); 
    neighbor_list_.CheckNeighborList(true);
    initialized_ = true;
}


void
ForceNeighborListCells::UpdateScheme() {
    // We check on each force calculation if we need to update
    // the neighbor list, so doesn't matter here.
    //neighbor_list_.CheckNeighborList();
}


void
ForceNeighborListCells::Interact() {
    if (!initialized_) {
        printf("ERROR: Finalized was not run for ForceNeighborListCells, exiting!\n");
        exit(1);
    }

    neighbor_list_.CheckNeighborList();


    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
        int tid;
        double **fr = new double*[3];
        double **tr = new double*[3];
        double *pr_energy;
        double *kmc_energy;

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
        kmc_energy = kmc_energy_ + (tid * nparticles_);
        for (int i = 0; i < nparticles_; ++i) {
            pr_energy[i] = 0.0;
            kmc_energy[i] = 0.0;
        }

        assert(nparticles_ == simples_->size());
        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Iterate over our neighbors
            for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
                int jdx = nldx->idx_;
                auto part1 = (*simples_)[idx];
                auto part2 = (*simples_)[jdx];

                // Do the interaction
                InteractParticlesMP(part1, part2, fr, tr, pr_energy, kmc_energy);
            }
        } // pragma omp for schedule(runtime) nowait

        // Reduce once all threads have finished
        #ifdef ENABLE_OPENMP
        #pragma omp barrier
        #endif
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
                kmc_energy_[jj] += kmc_energy_[offs + jj];
            }
        }

        delete[] fr;
        delete[] tr;
    } // pragma omp parallel

    ReduceParticlesMP();
}


// print specifics
void
ForceNeighborListCells::printSpecifics() {
    neighbor_list_.print();
}


// dump gory details
void
ForceNeighborListCells::dump() {
    ForceBase::dump();
    neighbor_list_.dump();
}
