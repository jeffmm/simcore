// Implementation for brute force calculation

#include "force_brute.h"

#include <cassert>

void
ForceBrute::InitMP() {
    // Nothing to do!
}


// Update the underlying scheme
// in this case, nothing!
void
ForceBrute::UpdateScheme() {
    // Nothing to do!
}


// Finalize everything, here it is fine, but have to set 
// the boolean flag
void
ForceBrute::Finalize() {
    initialized_ = true;
}


// Interactions between particles yay!
void
ForceBrute::Interact() {

    if (!initialized_) {
        printf("ERROR: Finalized was not run for ForceBrute, exiting!\n");
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
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
                auto part1 = simples_[idx];
                auto part2 = simples_[jdx];

                // Do the interaction itself from ForceBase
                InteractParticlesMP(part1, part2, fr, tr, pr_energy);
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
            }
        }

        delete[] fr;
        delete[] tr;
    } // pragma omp parallel

    ReduceParticlesMP();
}


void
ForceBrute::printSpecifics() {
    printf("\tNo specifics for brute force\n");
}
