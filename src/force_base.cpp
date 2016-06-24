// Implementation for the force_base

#include "force_base.h"

// Init the force
void
ForceBase::Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin) {
    space_ = pSpace;
    species_ = pSpecies;
    ndim_ = space_->n_dim;
    nperiodic_ = space_->n_periodic;
    skin_ = pSkin;
    for (int i = 0; i < ndim_; ++i) {
        box_[i] = space_->unit_cell[i][i];
    }

    #ifdef ENABLE_OPENMP
    #pragma omp parallel
    {
        if (0 == omp_get_thread_num()) {
            nthreads_ = omp_get_num_threads();
            //printf("Running OpenMP using %d threads\n", nthreads_);
        }
    }
    #else
    nthreads_ = 1;
    #endif
}


// XXX: CJE just copy the current scheme for adding potentials
// for now
void
ForceBase::InitPotentials(PotentialManager *pPotentials) {
    potentials_ = pPotentials;

    // Determine the maximum rcut
    max_rcut_ = potentials_->GetMaxRCut();
}


// Load the simple particles into the master vector
// Calculate the force superarray size
void
ForceBase::LoadSimples() {
   
    // Load everything from the species bases
    // Also remember that the OIDs may not be in order, so account for that
    // via a mapping oid <-> position
    nsys_ = species_->size();
    simples_.clear();
    for (auto it = species_->begin(); it != species_->end(); ++it) {
        std::vector<Simple*> sim_vec = (*it)->GetSimples();
        simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
    }
    nparticles_ = (int)simples_.size();

    // Ugh, figure out the OID stuff
    oid_position_map_.clear();
    for (int i = 0; i < nparticles_; ++i) {
        auto part = simples_[i];
        int oid = part->GetOID();
        // i is position in force superarray, oid is the actual oid
        oid_position_map_[oid] = i;
    }

    // Create the force and potential energy superarrays
    frc_ = new double[nthreads_*3*nparticles_];
    trqc_ = new double[nthreads_*3*nparticles_];
    prc_energy_ = new double[nthreads_*nparticles_];
}


// Main interaction routine for particles (once they have been determined by
// the force substructure
void
ForceBase::InteractParticlesMP(Simple *part1, Simple* part2, double **fr, double **tr, double *pr_energy) {
    // We are assuming the force/torque/energy superarrays are already set
    // XXX: JMM We should allow certain self-interactions
    // Exclude composite object interactions for now
    if (part1->GetCID() == part2->GetCID()) return;

    // Calculate the potential here
    PotentialBase *pot = potentials_->GetPotential(part1->GetSID(), part2->GetSID());
    if (pot == nullptr) return;
    // Minimum distance here@@@@!!!!
    // XXX: CJE ewwwwwww, more elegant way?
    interactionmindist idm;
    MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);
    if (idm.dr_mag2 > pot->GetRCut2()) return;

    // Obtain the mapping between particle oid and position in the force superarray
    auto oid1x = oid_position_map_[part1->GetOID()];
    auto oid2x = oid_position_map_[part2->GetOID()];
    #ifdef DEBUG
    if (debug_trace)
        printf("\tInteracting[%d:%d] (dr2:%2.2f)\n", oid1x, oid2x, idm.dr_mag2);
    #endif

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


// Reduce the particles back to their main versions
void
ForceBase::ReduceParticlesMP() {
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


// Print information about forces
void
ForceBase::print() {
    printf("********\n");
    printf("%s ->\n", name_.c_str());
    printf("\t{nthreads:%d}, {ndim:%d}, {nperiodic:%d}, {n:%d}\n", nthreads_, ndim_, nperiodic_, nparticles_);
    printf("\t{box:%2.2f}, {max_rcut:%2.2f}, {skin:%2.2f}\n", box_[0], max_rcut_, skin_);
    printSpecifics();
}


// Dump all of the information (very icky)
void
ForceBase::dump() {
    #ifdef DEBUG
    for (int i = 0; i < (int)simples_.size(); ++i) {
        auto part = simples_[i];
        auto oid = part->GetOID();
        printf("\to(%d) = ", oid);
        printf("x{%2.2f, %2.2f}, ", part->GetPosition()[0], part->GetPosition()[1]);
        printf("f{%2.2f, %2.2f}, ", part->GetForce()[0], part->GetForce()[1]);
        printf("u{%2.2f}, p{%2.2f}\n", part->GetKineticEnergy(), part->GetPotentialEnergy());
    }
    #endif
}


