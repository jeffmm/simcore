// Implementation for force microcells

#include <cassert>

#include "force_microcell.h"

// Overridden init method to call initmp
void
ForceMicrocell::Init(space_struct* pSpace, double pSkin) {
    // Override this to call base class, then initmp
    ForceBase::Init(pSpace, pSkin);
    InitMP();
}


// Overridden load simples method (needed for finalize)
void
ForceMicrocell::LoadSimples(std::vector<SpeciesBase*> pSpecies) {
    // Run the base class version
    ForceBase::LoadSimples(pSpecies);

    // Load the flat simples into the cell list
    microcell_list_.LoadFlatSimples(simples_);
}

void
ForceMicrocell::InitMP() {
    microcell_list_.Init(space_, skin_);
}


void
ForceMicrocell::Finalize() {
    // Should be ready to finalize if we have loaded and init'ed
    microcell_list_.CreateSubstructure(max_rcut_); 
    UpdateScheme();
    initialized_ = true;
}


void
ForceMicrocell::UpdateScheme() {
    microcell_list_.UpdateCellList();
}


void
ForceMicrocell::Interact() {
    if (!initialized_) {
        printf("ERROR: Finalized was not run for ForceMicrocells, exiting!\n");
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
        int ncells = microcell_list_.ncells();
        // Check within my own cell
        for (int cidx = 0; cidx < ncells; cidx += nthreads_) {
            // set the index
            int cjdx = cidx + tid;
            if (cjdx >= ncells) break;

            // Get the actual cell
            auto c1 = microcell_list_[cjdx];
            // Loop over particles in said cell
            for (int pidx1 = 0; pidx1 < c1->nparticles_ - 1; ++pidx1) {
                int ii = c1->idxlist_[pidx1]; // note,the location, not the oidx
                auto part1 = simples_[ii];

                // Get my interacting partner
                for (int pidx2 = pidx1 + 1; pidx2 < c1->nparticles_; ++pidx2) {
                    int jj = c1->idxlist_[pidx2];
                    auto part2 = simples_[jj];

                    // Do the interaction
                    InteractParticlesMP(part1, part2, fr, tr, pr_energy);
                }
            }
        } // cell loop


        // Interactions across different cells
        int npairs = microcell_list_.npairs();
        for (int pairidx = 0; pairidx < npairs; pairidx += nthreads_) {
            int pairjdx = pairidx + tid;
            if (pairjdx >= npairs) break;
            auto cell1 = microcell_list_[microcell_list_.plist(2*pairjdx  )];
            auto cell2 = microcell_list_[microcell_list_.plist(2*pairjdx+1)];

            for (int pidx1 = 0; pidx1 < cell1->nparticles_; ++pidx1) {
                int ii = cell1->idxlist_[pidx1];
                auto part1 = simples_[ii];

                for (int pidx2 = 0; pidx2 < cell2->nparticles_; ++pidx2) {
                    int jj = cell2->idxlist_[pidx2];
                    auto part2 = simples_[jj];

                    // Do the interaction
                    InteractParticlesMP(part1, part2, fr, tr, pr_energy);
                }
            }
        } // pair list

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


// print specific information
void
ForceMicrocell::printSpecifics() {
    microcell_list_.print();
}


// dump gory details
void
ForceMicrocell::dump() {
    ForceBase::dump();
    microcell_list_.dump();
}
