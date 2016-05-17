#include "system.h"

// Constructor
SystemArch::SystemArch(int nparticles, double pBox, double pDt) {
    next_pid_ = 0;

#if defined(_OPENMP)
#pragma omp parallel
    {
        if(0 == omp_get_thread_num()) {
            nthreads_ = omp_get_num_threads();
            printf("Runing OpenMP using %d threads\n", nthreads_);
        }
    }
#else
    nthreads_ = 1;
#endif
    nparticles_ = nparticles;
    box_ = pBox;
    ukin_ = 0.0;
    temperature_ = 0.0;
    dt_ = pDt;

    // Set up the main particle list
    particles_.clear();
    particles_.resize(nparticles_);

    // initialize the force superarray
    frc_.clear();
    frc_.resize(nthreads_ * 3 * nparticles_);
}


// Destructor
SystemArch::~SystemArch() {
    for(auto part : particles_) {
        delete(part);
    }
    for(auto spec : species_) {
        delete(spec.second);
    }
}


// Add a new species after construction
void
SystemArch::addSpecies(int sid, BaseSpecies* new_species) {
    species_[sid] = new_species;
    nsys_ = (int)species_.size();
}


// Add a new potential for pid1 and 2
void 
SystemArch::addPotential(int pid1, int pid2, PotentialBase* newPot) {
    potential_manager_.addPotential(pid1, pid2, newPot);
}


// Add a particle to species index idx
int
SystemArch::addParticle(int idx) {
    auto currentSpecies = getSpecies(idx);
    auto sid = currentSpecies->getSid();
    auto p  = currentSpecies->addParticle(sid, next_pid_);
    particles_[next_pid_] = p;
    next_pid_++;
    return next_pid_ - 1;
}


// Get the particular species in this slot
BaseSpecies*
SystemArch::getSpecies(int idx) {
    return species_[idx];
}


// Get the potential associted with pid1, 2
PotentialBase*
SystemArch::getPotential(int pid1, int pid2) {
    return potential_manager_.getPotential(pid1, pid2);
}


// Get the particle at index i
particle*
SystemArch::getParticle(int idx) {
    return particles_[idx];
}

// Flatten the particles into a single list
// Needed for cells...
void
SystemArch::flattenParticles() {
    // Get the total number of particles
    int nparticles_loc = 0;
    for(auto& sys : species_) {
        nparticles_loc += sys.second->getNParticles();
    }
    assert(nparticles_loc == nparticles_);
    
    // Organize the particles by sid
    std::sort(particles_.begin(), particles_.end());
}


// Generate the master cell list
void
SystemArch::generateCellList() {
    // Determine the total number of particles
    // And the box size
    // And the largest cutoff radius
    std::cout << "System generating cell list\n";
    double max_rcut = 0;

    for (int i = 0; i < nsys_; ++i) {
        auto currentSpecies = getSpecies(i);
        max_rcut = std::max(max_rcut, currentSpecies->getRcut());
    }

    // Flatten the particles
    flattenParticles();

    // Create the cell list
    double mbox[3] = {box_, box_, box_};
    cell_list_.CreateCellList(nparticles_, max_rcut, mbox);
    cell_list_.UpdateCellList(&particles_);
    cell_list_.CheckCellList();
}


// Update the cell list
void
SystemArch::updateCellList() {
    cell_list_.UpdateCellList(&particles_);
}


// Calculate the potential between two particles
void
SystemArch::calcPotential(int psid1, int psid2, double* x, double* y, double* fpote) {
    potential_manager_.calculatePotential(psid1, psid2, x, y, fpote);
}


// Calculate the kinetic energy
std::pair<double, double>
SystemArch::ukin() {
    std::pair<double, double> ukin = std::make_pair(0.0, 0.0);
    for (auto& sys : species_) {
        auto uk_temp = sys.second->Ukin(&particles_);
        ukin.first += uk_temp.first;
        ukin.second += uk_temp.second;
    }
    ukin_ = ukin.first;
    temperature_ = ukin.second;
    return ukin;
}


// MP Force calculation routine
void
SystemArch::forceMP() {
    double epot = 0.0;

#if defined(_OPENMP)
#pragma omp parallel reduction(+:epot)
#endif
    {
        int tid;
        double *fx, *fy, *fz;
        double f_epot[4];

#if defined(_OPENMP)
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif

        // Set up the pointers to the force superarray
        fx = frc_.data() + (3*tid*nparticles_);
        buffmd::azzero(fx, 3*nparticles_);
        fy = frc_.data() + ((3*tid+1)*nparticles_);
        fz = frc_.data() + ((3*tid+2)*nparticles_);

        // Check within my own cell
        int ncells = cell_list_.ncells();
        for (int cidx = 0; cidx < ncells; cidx += nthreads_) {
            // set the index
            int cjdx = cidx + tid;
            if (cjdx >= ncells) break;

            // Get the actual cell
            auto c1 = cell_list_[cjdx];
            // Loop over particles in said cell
            for (int pidx1 = 0; pidx1 < c1->nparticles_ - 1; ++pidx1) {
                int ii = c1->idxlist_[pidx1];
                auto part1 = particles_[ii];

                // Get my interacting partner
                for(int pidx2 = pidx1 + 1; pidx2 < c1->nparticles_; ++pidx2) {
                    int jj = c1->idxlist_[pidx2];
                    auto part2 = particles_[jj];

                    calcPotential(part1->sid, part2->sid, part1->x, part2->x, f_epot);
                    epot += f_epot[3];
                    fx[ii] += f_epot[0];
                    fy[ii] += f_epot[1];
                    fz[ii] += f_epot[2];
                    fx[jj] -= f_epot[0];
                    fy[jj] -= f_epot[1];
                    fz[jj] -= f_epot[2];
                } // check interaction partner particles
            } // check the actual particles
        }  // Cell loop

        // Interactions across different cells
        int npairs = cell_list_.npairs();
        for (int pairidx = 0; pairidx < npairs; pairidx += nthreads_) {
            int pairjdx = pairidx + tid;
            if (pairjdx >= npairs) break;
            auto cell1 = cell_list_[cell_list_.plist(2*pairjdx  )];
            auto cell2 = cell_list_[cell_list_.plist(2*pairjdx+1)];

            for (int pidx1 = 0; pidx1 < cell1->nparticles_; ++pidx1) {
                int ii = cell1->idxlist_[pidx1];
                auto part1 = particles_[ii];

                for (int pidx2 = 0; pidx2 < cell2->nparticles_; ++pidx2) {
                    int jj = cell2->idxlist_[pidx2];
                    auto part2 = particles_[jj];

                    calcPotential(part1->sid, part2->sid, part1->x, part2->x, f_epot);
                    epot += f_epot[3];
                    fx[ii] += f_epot[0];
                    fy[ii] += f_epot[1];
                    fz[ii] += f_epot[2];
                    fx[jj] -= f_epot[0];
                    fy[jj] -= f_epot[1];
                    fz[jj] -= f_epot[2];
                } // Second particle
            } // First particle
        } // Pairs loops

        // reduce once all threads have finished
#if defined(_OPENMP)
#pragma omp barrier
#endif
        int i = 1 + (3 * nparticles_ / nthreads_);
        int fromidx = tid * i;
        int toidx = fromidx + i;
        if (toidx > 3*nparticles_) toidx = 3*nparticles_;

        // Reduce the forces
        for (i = 1; i < nthreads_; ++i) {
            int offs;

            offs = 3*i*nparticles_;

            for (int j = fromidx; j < toidx; ++j) {
                frc_[j] += frc_[offs+j];
            }
        }
    } // omp parallel reduction epot

    // Recombine into the particles
    for (int i = 0; i < nparticles_; ++i) {
        auto part = particles_[i];
        part->f[0] = frc_[i];
        part->f[1] = frc_[nparticles_ + i];
        part->f[2] = frc_[2*nparticles_ + i];
    }
    upot_ = epot;
}


// Position and velocity update
void
SystemArch::velverlet() {
    // Update initial step and velocities
    for (int idx = 0; idx < nparticles_; ++idx) {
        auto p = particles_[idx];
        auto pmeff = species_[p->sid]->getMeff();
        double dtmf = 0.5 * dt_ / pmeff;
        for (int i = 0; i < 3; ++i) {
            p->v[i] += dtmf * p->f[i];
            p->x[i] += dt_ * p->v[i];
        }
    }

    // Compute energies
    forceMP();

    // Update another half step
    for (int idx = 0; idx < nparticles_; ++idx) {
        auto p = particles_[idx];
        auto pmeff = species_[p->sid]->getMeff();
        double dtmf = 0.5 * dt_ / pmeff;
        for (int i = 0; i < 3; ++i) {
            p->v[i] += dtmf * p->f[i];
        }
    }

}


// output the simluation information
void
SystemArch::output(FILE* erg, FILE* traj, int nfi) {
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", nfi, temperature_, ukin_, upot_, ukin_+upot_);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", nfi, temperature_, ukin_, upot_, ukin_+upot_);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", nparticles_, nfi, ukin_+upot_);
    // Loop over base particles
    for (int i=0; i < nparticles_; ++i) {
        auto p = particles_[i];
        fprintf(traj, "%s  %20.8f %20.8f %20.8f\n", p->name.c_str(), p->x[0], p->x[1], p->x[2]);
    }
}


// Get the total number of particles
int
SystemArch::nParticles() {
    return nparticles_;
}


// Dump EVERYTHING!!
void
SystemArch::dump() {
    std::cout << "********\n";
    std::cout << "System dump: \n";
    printf("{nparticles: %d}, {box: %f}\n", nparticles_, box_);
    for(auto sys : species_) {
        std::cout << "Species: ";
        sys.second->dump();
    }
}


// Dump the potentials information
void
SystemArch::dumpPotentials() {
    std::cout << "********\n";
    std::cout << "Potentials dump: \n";
    potential_manager_.print();
}


// Check the consistency of everything
void
SystemArch::checkConsistency() {
    std::cout << "********\n";
    std::cout << "Consistency Check: \n";
    for(auto sys : species_) {
        sys.second->checkParticles();
    }
}


