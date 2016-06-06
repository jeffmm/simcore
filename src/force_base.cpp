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

    #ifdef ENABLE_OPNEMP
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
        for (auto jt=pot_vec.begin(); jt!=pot_vec.end(); ++jt)
            potentials_.AddPotential(jt->first.first,jt->first.second,jt->second);
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
        auto oid = part->GetOID();
        // i is position in force superarray, oid is the actual oid
        oid_position_map_[oid] = i;
    }

    // Resize the force and potential energy superarrays
    frc_.clear();
    frc_.resize(nthreads_ * ndim_ * nparticles_);
    prc_energy_.clear();
    prc_energy_.resize(nthreads_ * nparticles_);
}


// Find the minimum distance beween two particles
void ForceBase::MinimumDistance(Simple* o1, Simple* o2, interactionmindist& imd) {
    double const * const r1 = o1->GetPosition();
    double const * const s1 = o1->GetScaledPosition();
    double const * const u1 = o1->GetOrientation();
    double const * const r2 = o2->GetPosition();
    double const * const s2 = o2->GetScaledPosition();
    double const * const u2 = o2->GetOrientation();
    double const l1 = o1->GetLength();
    double const l2 = o2->GetLength();
    double const d1 = o1->GetDiameter();
    double const d2 = o2->GetDiameter();
    /* TODO: Think about how best to do this for general shapes, like 2d
       polygons that can represent the local surface of more complex 3d
       shapes. Perhaps assume all local surface to be triangular polygons.*/
    imd.dr_mag2 = 0;
    std::fill(imd.dr, imd.dr+3, 0.0);
    std::fill(imd.contact1, imd.contact1+3, 0.0);
    std::fill(imd.contact2, imd.contact2+3, 0.0);
    imd.buffer_mag = 0.5*(d1+d2);
    imd.buffer_mag2 = imd.buffer_mag*imd.buffer_mag;
    if (l1 == 0 && l2 == 0)
        min_distance_point_point(ndim_, nperiodic_, space_->unit_cell, 
                                 r1, s1, r2, s2, imd.dr, &imd.dr_mag2);
    else if (l1 == 0 && l2 > 0) 
        min_distance_sphere_sphero(ndim_, nperiodic_, space_->unit_cell,
                                   r1, s1, r2, s2, u2, l2,
                                   imd.dr, &imd.dr_mag2, imd.contact2);
    else if (l1 > 0 && l2 == 0) 
        min_distance_sphere_sphero(ndim_, nperiodic_, space_->unit_cell,
                                   r2, s2, r1, s1, u1, l1,
                                   imd.dr, &imd.dr_mag2, imd.contact1);
    else if (l1 > 0 && l2 > 0)
        min_distance_sphero(ndim_, nperiodic_, space_->unit_cell,
                            r1, s1, u1, l1, r2, s2, u2, l2,
                            imd.dr, &imd.dr_mag2, imd.contact1, imd.contact2);
    imd.dr_mag = sqrt(imd.dr_mag2);
}
