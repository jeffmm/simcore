// Implementation for interaction engine

#include "interaction_engine.h"

void InteractionEngine::Init(space_struct *pSpace, ParticleTracking *pTracking, double pSkin) {
  space_ = pSpace;
  tracking_ = pTracking;
  skin_= pSkin;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  for (int i = 0; i < ndim_; ++i) {
    box_[i] = space_->unit_cell[i][i];
  }
  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  {
    if (0 == omp_get_thread_num()) {
      nthreads_ = omp_get_num_threads();
      printf("Running with %d OpenMP threads\n", nthreads_);
    }
  }
  #else
  nthreads_ = 1;
  #endif
}

// Initialize the potentials, get the maximum rcut value
void InteractionEngine::InitPotentials(PotentialManager *pPotentials) {
  potentials_ = pPotentials;
  max_rcut_ = potentials_->GetMaxRCut();
}

// Initialize the stuff needed for MP (superarrays, et al)
void InteractionEngine::InitMP() {
  nsimples_ = tracking_->GetNSimples();
  simples_ = tracking_->GetSimples();

  // Clear out things from before (if they exist)
  if (frc_ != nullptr) {
    delete[] frc_;
  }
  if (trqc_ != nullptr) {
    delete[] trqc_;
  }
  if (prc_energy_ != nullptr) {
    delete[] prc_energy_;
  }
  if (kmc_energy_ != nullptr) {
    delete[] kmc_energy_;
  }

  // Ugh, figure out the OID stuff
  oid_position_map_.clear();
  for (int i = 0; i < nsimples_; ++i) {
    auto part = (*simples_)[i];
    int oid = part->GetOID();
    oid_position_map_[oid] = i;
  }

  // Create the force and potential energy superarrays
  frc_ = new double[nthreads_*3*nsimples_];
  trqc_ = new double[nthreads_*3*nsimples_];
  prc_energy_ = new double[nthreads_*nsimples_];
  kmc_energy_ = new double[nthreads_*nsimples_];
}

// Actual interaction routine
// Always uses neighbor list
void InteractionEngine::Interact() {

  // XXX Check reinit?

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    int tid;
    double **fr = new double*[3];
    double **tr = new double*[3];
    double *pr_energy;
    double *kmc_energy;

    auto neighbors = tracking_->GetNeighbors();

    #ifdef ENABLE_OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    fmmph::InitMPRegion(tid, nsimples_, &frc_, &trqc_, &prc_energy_, &kmc_energy_,
        &fr, &tr, &pr_energy, &kmc_energy);

    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int idx = 0; idx < nsimples_; ++idx) {
      // Iterate over our neighbors
      for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
        int jdx = nldx->idx_;
        auto part1 = (*simples_)[idx];
        auto part2 = (*simples_)[jdx];

        // Do the interaction
        InteractParticlesMP(&(*nldx), part1, part2, fr, tr, pr_energy, kmc_energy);
      }
    } // pragma omp for schedule(runtime) nowait

    // Reduce once all threads have finished
    #ifdef ENABLE_OPENMP
    #pragma omp barrier
    #endif
    fmmph::ReduceMPRegion(tid, nsimples_, nthreads_, &frc_, &trqc_, &prc_energy_, &kmc_energy_);

    delete[] fr;
    delete[] tr;
  } // pragma omp parallel
  ReduceParticlesMP();
}

// Main interaction routine for particles 
void InteractionEngine::InteractParticlesMP(neighbor_t *neighbor, Simple *part1, Simple* part2, double **fr, double **tr, double *pr_energy, double *kmc_energy) {
  // We are assuming the force/torque/energy superarrays are already set
  // Exclude composite object interactions
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

  // Switch on whether or not we are doing a mc or force calculation
  if (pot->IsKMC()) {
    // Figure out which particle is the kmc particle
    SID kmc_target = pot->GetKMCTarget();
    if (kmc_target == part1->GetSID()) {
      kmc_energy[oid1x] += fepot[ndim_];
    }
    if (kmc_target == part2->GetSID()) {
      kmc_energy[oid2x] += fepot[ndim_];
    }
  } else {
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
}

// Reduce the particles back to their main versions
void InteractionEngine::ReduceParticlesMP() {
  for (int i = 0; i < nsimples_; ++i) {
    auto part = (*simples_)[i];
    int oidx = oid_position_map_[part->GetOID()];
    double subforce[3] = {0.0, 0.0, 0.0};
    double subtorque[3] = {0.0, 0.0, 0.0};
    for (int idim = 0; idim < ndim_; ++idim) {
      subforce[idim] = frc_[idim*nsimples_+oidx];
      subtorque[idim] = trqc_[idim*nsimples_+oidx]; 
    }
    part->AddForceTorqueEnergyKMC(subforce, subtorque, prc_energy_[oidx], kmc_energy_[oidx]);
  }
}

// Print out information
void InteractionEngine::Print() {
  printf("********\n");
  printf("InteractionEngine ->\n");
  printf("\t{nthreads:%d}, {ndim:%d}, {nperiodic:%d}, {n:%d}\n", nthreads_, ndim_, nperiodic_, nsimples_);
  printf("\t{box:%2.2f}, {max_rcut:%2.2f}, {skin:%2.2f}\n", box_[0], max_rcut_, skin_);
}

// Dump all the glorious information
void InteractionEngine::Dump() {
  #ifdef DEBUG
  if (debug_trace) {
    for (int i = 0; i < nsimples_; ++i) {
      auto part = (*simples_)[i];
      auto oid = part->GetOID();
      printf("\to(%d) = ", oid);
      printf("x{%2.2f, %2.2f}, ", part->GetPosition()[0], part->GetPosition()[1]);
      printf("f{%2.2f, %2.2f}, ", part->GetForce()[0], part->GetForce()[1]);
      printf("u{%2.2f}, p{%2.2f}\n", part->GetKineticEnergy(), part->GetPotentialEnergy());
    }
  }
  #endif
}
