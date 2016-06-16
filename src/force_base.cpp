// Implementation for the force_base

#include "force_base.h"

// Init the force
void
ForceBase::Init(space_struct *pSpace, double pSkin) {
    space_ = pSpace;
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
            printf("Running OpenMP using %d threads\n", nthreads_);
        }
    }
    #else
    nthreads_ = 1;
    #endif
}


// XXX: CJE just copy the current scheme for adding potentials
// for now
void
ForceBase::InitPotentials(std::vector<SpeciesBase*> pSpecies) {
    for (auto it = pSpecies.begin(); it != pSpecies.end(); ++it) {
        auto pot_vec = (*it)->GetPotentials();
        for (auto jt=pot_vec.begin(); jt!=pot_vec.end(); ++jt) {
            potentials_.AddPotential(jt->first.first,jt->first.second,jt->second);
            max_rcut_ = std::max(max_rcut_, jt->second->GetRCut());
        }
    }
    
    potentials_.Print();
}


// Load the simple particles into the master vector
// Calculate the force superarray size
void
ForceBase::LoadSimples(std::vector<SpeciesBase*> pSpecies) {
   
    // Load everything from the species bases
    // Also remember that the OIDs may not be in order, so account for that
    // via a mapping oid <-> position
    nsys_ = pSpecies.size();
    simples_.clear();
    for (auto it = pSpecies.begin(); it != pSpecies.end(); ++it) {
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

