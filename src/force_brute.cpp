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


// Interactions between particles yay!
void
ForceBrute::Interact() {
 
    printf("Brute force interact starting parallel section\n");
    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    #endif
    {
        int tid;
        double **fr;
        fr = (double**)malloc(ndim_ * sizeof(double*));
        double *pr_energy;

        #ifdef ENABLE_OPENMP
        tid = omp_get_thread_num();
        #else
        tid = 0;
        #endif

        // Set up the pointers to the force superarray
        for (int i = 0; i < ndim_; ++i) {
            fr[i] = frc_.data() + ((ndim_*tid+i)*nparticles_);
        }
        std::fill(&fr[0][0], &fr[0][0] + ndim_*nparticles_, 0.0);
        pr_energy = prc_energy_.data() + (tid * nparticles_);
        std::fill(pr_energy, pr_energy + nparticles_, 0.0);

        printf("\t(%d) -> n(%d), ndim(%d), starting loop\n", tid, nparticles_, ndim_);

        assert(nparticles_ == simples_.size());
        #ifdef ENABLE_OPENMP
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
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
                MinimumDistance(part1, part2, idm);
                if (idm.dr_mag2 > pot->GetRCut2()) continue;

                // Obtain the mapping between particle oid and position in the force superarray
                auto oid1x = oid_position_map_[part1->GetOID()];
                auto oid2x = oid_position_map_[part2->GetOID()];

                // Fire off the potential calculation
                double fepot[4];
                pot->CalcPotential(idm.dr, idm.dr_mag, idm.buffer_mag, fepot);
                pr_energy[oid1x] += fepot[ndim_];
                pr_energy[oid2x] += fepot[ndim_];
                for (int i = 0; i < ndim_; ++i) {
                    fr[i][oid1x] += fepot[i];
                    fr[i][oid2x] -= fepot[i];
                }
            }
        } // pragma omp for schedule(runtime) nowait

        // Reduce once all threads have finished
        #ifdef ENABLE_OPENMP
        #pragma omp barrier
        #endif
        printf("\t(%d) -> done with loop, reducing!\n", tid);
        // Everything is doubled to do energy
        // XXX: CJE maybe fix this later?
        int i = 1 + (ndim_ * nparticles_ / nthreads_);
        int ii = 1 + (nparticles_ / nthreads_);

        int fromidx = tid * i;
        int fromidx2 = tid * ii;

        int toidx = fromidx + i;
        int toidx2 = fromidx2 + ii;

        if (toidx > ndim_*nparticles_) toidx = ndim_*nparticles_;
        if (toidx2 > nparticles_) toidx2 = nparticles_;

        // Reduce the forces
        for (i = 1; i < nthreads_; ++i) {
            int offs;
            offs = ndim_ * i * nparticles_;

            for (int j = fromidx; j < toidx; ++j) {
                frc_[j] += frc_[offs+j];
            }
        }

        // Reduce the energies
        for (ii = 1; ii < nthreads_; ++ii) {
            int offs;
            offs = i * nparticles_;

            for (int jj = fromidx2; jj < toidx2; ++jj) {
                prc_energy_[jj] += prc_energy_[offs + jj];
            }
        }

        delete(fr);
    } // pragma omp parallel

    printf("Brute force interact done with parallel section\n");

    // XXX: CJE this is where we tell the particles what to do
    for (int i = 0; i < nparticles_; ++i) {
        auto part = simples_[i];
        auto oidx = oid_position_map_[part->GetOID()];
        double subforce[3] = {frc_[oidx], frc_[nparticles_+oidx], frc_[2*nparticles_+oidx]};
        // XXX: CJE Replace with actual torque
        double subtorque[3] = {frc_[oidx], frc_[nparticles_+oidx], frc_[2*nparticles_+oidx]};
        part->AddForceTorqueEnergy(subforce, subtorque, prc_energy_[oidx]);
    }

    printf("Brute force interact returning\n");
}
