// Implementation for force adjacent cells

#include <cassert>

#include "force_cell.h"

// Overridden init method to call initmp
void
ForceCell::Init(space_struct* pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin) {
    // Override this to call base class, then initmp
    ForceBase::Init(pSpace, pSpecies, pSkin);
    InitMP();
}


// Overridden load simples method (needed for finalize)
void
ForceCell::LoadSimples() {
    // Run the base class version
    ForceBase::LoadSimples();

    // Load the flat simples into the cell list
    cell_list_.LoadFlatSimples(simples_);
}

void
ForceCell::InitMP() {
    cell_list_.Init(space_, species_, skin_);
}


void
ForceCell::Finalize() {
    // Should be ready to finalize if we have loaded and init'ed
    cell_list_.CreateSubstructure(max_rcut_); 
    UpdateScheme();
    initialized_ = true;
}


void
ForceCell::UpdateScheme() {
    cell_list_.UpdateCellList();
}


void
ForceCell::Interact() {
    if (!initialized_) {
        printf("ERROR: Finalized was not run for ForceCells, exiting!\n");
        exit(1);
    }

    // Do the interactions!
    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
        int tid;
        double **fr = new double*[3];
        double **tr = new double*[3];
        double *pr_energy;

        std::vector<int>* pid_to_cid = cell_list_.pidtocid();

        #ifdef ENABLE_OPENMP
        tid = omp_get_thread_num();
        #else
        tid = 0;
        #endif

        // Set up pointers to the superarrays
        for (int i = 0; i < 3; ++i) {
            fr[i] = frc_ + ((3*tid+i)*nparticles_);
            tr[i] = trqc_ + ((3*tid+i)*nparticles_);
        }

        for (int i = 0; i < 3*nparticles_; ++i) {
            (*fr)[i] = 0.0;
            (*tr)[i] = 0.0;
        }

        pr_energy = prc_energy_ + (tid*nparticles_);
        for (int i = 0; i < nparticles_; ++i) {
            pr_energy[i] = 0.0;
        }

        assert(nparticles_ == simples_.size());
        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Get our cell
            int cidx = (*pid_to_cid)[idx];
            auto cell1 = cell_list_[cidx];
            // Loop over all other cells (including this one)
            int nadj = cell_list_.nadj();
            for (int cjdx = 0; cjdx < nadj; ++cjdx) {
                auto cell2 = cell_list_[cell1->adj_cell_ids_[cjdx]];
                // Loop over the second cell's particles
                for (int jdx = 0; jdx < cell2->nparticles_; ++jdx) {
                    int jjdx = cell2->idxlist_[jdx];
                    // ONLY DO THE CALCULATION IF THE OTHER PID IS HIGHER!!!
                    if (jjdx > idx) {
                        auto part1 = simples_[idx];
                        auto part2 = simples_[jjdx];

                        // Interact
                        InteractParticlesMP(part1, part2, fr, tr, pr_energy);
                    } // only do this calculation if the second pid is higher
                } // cell2 particle list
            } // cells adjacent  and equal to us
        } // pragma omp for schedule(runtime) nowait

        // Reduce once all threads have finished
        #ifdef ENABLE_OPENMP
        #pragma omp barrier
        #endif
        // Everything doubled to do energy and torque
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

        delete[] fr;
        delete[] tr;
    } // pragma omp parallel

    ReduceParticlesMP();
}


// print specifics
void
ForceCell::printSpecifics() {
    cell_list_.print();
}


// dump gory details
void
ForceCell::dump() {
    ForceBase::dump();
    cell_list_.dump();
}
