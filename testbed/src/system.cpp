#include <cstdlib>

#include "system.h"


// Constructor
SystemArch::SystemArch(properties_t* pProperties) {
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
    system_properties_ = pProperties;
    nparticles_ = system_properties_->nparticles_;
    ndim_ = system_properties_->ndim_;
    memcpy(box_, system_properties_->box_, 3*sizeof(double));
    skin_ = system_properties_->skin_;
    // Make sure that we don't accidentally use the wrong skin for
    // the cell list (unless we're using it as part of a neighbor list
    if (system_properties_->scheme_ == FCELLS) {
        skin_ = 0.0;
    }
    dt_ = system_properties_->dt_;
    
    ukin_ = 0.0;
    temperature_ = 0.0;
    
    // Set up the main particle list
    particles_.clear();
    particles_.resize(nparticles_);
    
    // initialize the force superarray
    frc_.clear();
    frc_.resize(nthreads_ * ndim_ * nparticles_);
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


// Initialize the shared data structures, etc
void
SystemArch::initMP() {
    printf("********\n");
    printf("Initializing MP shared data structures\n");
    flattenParticles();
    switch(system_properties_->scheme_) {
        case BRUTEFORCE:
            break;
        case FCELLS:
            generateCellList();
            break;
        case FNEIGHBORS_ALLPAIRS:
            generateNeighborList();
            break;
        case FNEIGHBORS_CELL:
            generateNeighborListCell();
            break;
        default:
            fprintf(stderr,"Not a supported data structure type!\n");
            exit(1);
    }
}


// Generate the master cell list
void
SystemArch::generateCellList() {
    // Determine the total number of particles
    // And the box size
    // And the largest cutoff radius
    double max_rcut = 0;

    for (int i = 0; i < nsys_; ++i) {
        auto currentSpecies = getSpecies(i);
        max_rcut = std::max(max_rcut, currentSpecies->getRcut());
    }

    // Create the cell list
    cell_list_.CreateCellList(nparticles_, max_rcut, skin_, box_);
    cell_list_.UpdateCellList(&particles_);
    cell_list_.CheckCellList();
}


// Update the cell list
void
SystemArch::updateCellList() {
    cell_list_.UpdateCellList(&particles_);
}


// Generate the neighbor list
void
SystemArch::generateNeighborList() {
    double max_rcut = 0.0;
    
    for (int i = 0; i < nsys_; ++i) {
        auto current_species = getSpecies(i);
        max_rcut = std::max(max_rcut, current_species->getRcut());
    }
    
    neighbor_list_.CreateNeighborList(nparticles_, max_rcut, skin_, box_);
    neighbor_list_.UpdateNeighborList(&particles_);
    neighbor_list_.print();
}


// Generate the neighbor list (cells)
void
SystemArch::generateNeighborListCell() {
    double max_rcut = 0.0;
    
    for (int i = 0; i < nsys_; ++i) {
        auto current_species = getSpecies(i);
        max_rcut = std::max(max_rcut, current_species->getRcut());
    }
    
    neighbor_list_cell_.CreateNeighborList(nparticles_, max_rcut, skin_, box_);
    neighbor_list_cell_.SetNThreads(nthreads_);
    neighbor_list_cell_.UpdateNeighborList(&particles_);
    neighbor_list_cell_.print();
}


// Calculate the potential between two particles
void
SystemArch::calcPotential(int psid1, int psid2, double* x, double* y, double* fpote) {
    potential_manager_.calculatePotential(psid1, psid2, x, y, fpote);
}


// Calculate the kinetic energy and temperature
std::pair<double, double>
SystemArch::ukin() {
    ukin_ = 0.0;
    temperature_ = 0.0;
    // Get the kinetic energy of all particles
    // Then calculate the temperature
    for (auto& sys : species_) {
        ukin_ += sys.second->Ukin(&particles_);
    }
    temperature_ = 2.0 * ukin_ / (ndim_ * nparticles_ - ndim_)/kboltz;
    //temperature_ = 2.0 * ukin_ / (3.0 * nparticles_ - 3.0)/kboltz;
    return std::make_pair(ukin_, temperature_);
}


// Brute force calculation routine, for testing purposes ONLY!!!
void
SystemArch::forceBrute() {
    
    double epot = 0.0;
    
#if defined(_OPENMP)
#pragma omp parallel
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
        
#if defined(_OPENMP)
#pragma omp for reduction(+:epot) schedule(runtime) nowait
#endif
        for (int idx = 0; idx < nparticles_ - 1; ++idx) {
            for (int jdx = idx + 1; jdx < nparticles_; ++jdx) {
                auto part1 = particles_[idx];
                auto part2 = particles_[jdx];
                // Calculate the potential (takes care of cutoff)
                
                calcPotential(part1->sid, part2->sid, part1->x, part2->x, f_epot);
                epot += f_epot[3];
                fx[idx] += f_epot[0];
                fy[idx] += f_epot[1];
                fz[idx] += f_epot[2];
                fx[jdx] -= f_epot[0];
                fy[jdx] -= f_epot[1];
                fz[jdx] -= f_epot[2];
            }
        } // pragma omp for reduction(+:epot) schedule(runtime) nowait
        
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
    } // pragma omp parallel
    
    // Recombine into the particles
    for (int i = 0; i < nparticles_; ++i) {
        auto part = particles_[i];
        part->f[0] = frc_[i];
        part->f[1] = frc_[nparticles_ + i];
        part->f[2] = frc_[2*nparticles_ + i];
    }
    upot_ = epot;
}


// Force calculation routine (in general, will replace forceMP)
void
SystemArch::forceNeighAP() {
    
    double epot = 0.0;
    
    // check the neighbor list for updates!
    neighbor_list_.CheckNeighborList(&particles_);
    
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
        int tid;
        double *fx, *fy, *fz;
        double f_epot[4];
        auto neighbors = neighbor_list_.GetNeighbors();
        
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

#if defined(_OPENMP)
#pragma omp for reduction(+:epot) schedule(runtime) nowait
#endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Iterate over the entries in our neighbor list
            for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); nldx++) {
                int jdx = nldx->idx_;
                auto part1 = particles_[idx];
                auto part2 = particles_[jdx];
                // Calculate the potential (takes care of cutoff)
                
                calcPotential(part1->sid, part2->sid, part1->x, part2->x, f_epot);
                epot += f_epot[3];
                fx[idx] += f_epot[0];
                fy[idx] += f_epot[1];
                fz[idx] += f_epot[2];
                fx[jdx] -= f_epot[0];
                fy[jdx] -= f_epot[1];
                fz[jdx] -= f_epot[2];
            }
        } // pragma omp for reduction(+:epot) schedule(runtime) nowait
        
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
    } // pragma omp parallel
    
    // Recombine into the particles
    for (int i = 0; i < nparticles_; ++i) {
        auto part = particles_[i];
        part->f[0] = frc_[i];
        part->f[1] = frc_[nparticles_ + i];
        part->f[2] = frc_[2*nparticles_ + i];
    }
    upot_ = epot;
}


// Force calculation routine (in general, will replace forceMP)
void
SystemArch::forceNeighCell() {
    
    double epot = 0.0;
    
    // check the neighbor list for updates!
    neighbor_list_cell_.CheckNeighborList(&particles_);
    
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
        int tid;
        double *fx, *fy, *fz;
        double f_epot[4];
        auto neighbors = neighbor_list_cell_.GetNeighbors();
        
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
        
#pragma omp for reduction(+:epot) schedule(runtime) nowait
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Iterate over the entries in our neighbor list
            for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); nldx++) {
                int jdx = nldx->idx_;
                auto part1 = particles_[idx];
                auto part2 = particles_[jdx];
                // Calculate the potential (takes care of cutoff)
                
                calcPotential(part1->sid, part2->sid, part1->x, part2->x, f_epot);
                epot += f_epot[3];
                fx[idx] += f_epot[0];
                fy[idx] += f_epot[1];
                fz[idx] += f_epot[2];
                fx[jdx] -= f_epot[0];
                fy[jdx] -= f_epot[1];
                fz[jdx] -= f_epot[2];
            }
        } // pragma omp for reduction(+:epot) schedule(runtime) nowait
        
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
    } // pragma omp parallel
    
    // Recombine into the particles
    for (int i = 0; i < nparticles_; ++i) {
        auto part = particles_[i];
        part->f[0] = frc_[i];
        part->f[1] = frc_[nparticles_ + i];
        part->f[2] = frc_[2*nparticles_ + i];
    }
    upot_ = epot;
}


// MP Force calculation routine
void
SystemArch::forceCellsMP() {
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
            // Accumulator for the neighbor list update
            p->dr_tot[i] += dt_ * p->v[i];
        }
    }

    // Compute energies
    switch(system_properties_->scheme_) {
        case BRUTEFORCE:
            forceBrute();
            break;
        case FCELLS:
            forceCellsMP();
            break;
        case FNEIGHBORS_ALLPAIRS:
            forceNeighAP();
            break;
        case FNEIGHBORS_CELL:
            forceNeighCell();
            break;
        default:
            fprintf(stderr,"Something has gone horribly wrong!\n");
            exit(1);
    }

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


// Calculate end statistics from this run
void
SystemArch::statistics(int pNsteps) {
    int nupdates = neighbor_list_.GetNUpdates();
    printf("Neighbor List:\n");
    printf("\t%d/%d updates/steps\n", nupdates, pNsteps);
}


// Dump EVERYTHING!!
void
SystemArch::dump() {
    std::cout << "********\n";
    std::cout << "System dump: \n";
    printf("{nparticles: %d}, {box: (%f,%f,%f)}\n", nparticles_, box_[0], box_[1], box_[2]);
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


