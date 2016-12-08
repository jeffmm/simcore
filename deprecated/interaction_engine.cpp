// Implementation of new interaction engine

#include "interaction_engine.h"

#include <cassert>

void InteractionEngine::Init(space_struct *pSpace,
                               ParticleEngine *pTrackEngine,
                               std::vector<interaction_t> *pInteractions) {
  space_ = pSpace;
  ptrack_ = pTrackEngine;
  interactions_ = pInteractions;

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
      std::cout << "Running with " << nthreads_ << " OpenMP threads\n";
    }
  }
  #else
  nthreads_ = 1;
  #endif

  AttachParticleEngine();
}

void InteractionEngine::AttachParticleEngine() {
  // Attach to the particle engine
  simples_  = ptrack_->GetSimples();
  species_  = ptrack_->GetSpecies();
  oid_position_map_ = ptrack_->GetOIDPositionMap();
  anchors_  = ptrack_->GetAnchors();
}

void InteractionEngine::InitMP() {
  nsimples_ = (int)simples_->size();
  nspecies_ = (int)species_->size();

  for (int i = 0; i < nspecies_; ++i) {
    spec_ind_map_[(*species_)[i]->GetSID()] = i;
  }

  if (frc_ != nullptr) {
    delete[] frc_;
  }
  if (trqc_ != nullptr) {
    delete[] trqc_;
  }
  if (prc_energy_ != nullptr) {
    delete[] prc_energy_;
  }
  if (virial_ != nullptr) {
    delete[] virial_;
  }

  // Create the superarrays
  frc_  = new double[nthreads_*3*nsimples_];
  trqc_ = new double[nthreads_*3*nsimples_];
  prc_energy_ = new double[nthreads_*nsimples_];
  virial_     = new double[nthreads_*9*nspecies_];
}

void InteractionEngine::Interact() {
  // Particle engine should have already updated if necessary
  // Update the other information of note

  int new_nsimples  = (int)simples_->size();
  int new_nspecies  = (int)species_->size();
  // If we've had a change in particle number or something, then
  // we need to update the superarrays
  if ((new_nsimples != nsimples_) ||
      (new_nspecies != nspecies_)) {
    std::cout << "Interaction engine got a different number of simples, rebuilding\n";
    InitMP();
  }

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    int tid;
    double **fr = new double*[3];
    double **tr = new double*[3];
    double **virial = new double*[9];
    double *pr_energy;

    #ifdef ENABLE_OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    // Set up the connection to the superarrays
    ieh::InitMPRegion(tid,
                      nsimples_,
                      nspecies_,
                      &frc_,
                      &trqc_,
                      &prc_energy_,
                      &virial_,
                      &fr,
                      &tr,
                      &pr_energy,
                      &virial);

    // Loop over the interactions
    int nix = (int)interactions_->size();
    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int ix = 0; ix < nix; ++ix) {
      // Alias the interaction
      interaction_t *pix = &((*interactions_)[ix]);
      //std::cout << "[" << ix << "] -> {" << pix->idx_ << " -> " << pix->jdx_ << "}\n";
      switch(pix->type_) {
        case ptype::external:
          InteractParticlesExternalMP(&pix, fr, tr, pr_energy, virial);
          break;
        case ptype::kmc:
          InteractParticlesKMCMP(&pix, fr, tr, pr_energy, virial);
          break;
        case ptype::internal:
          InteractParticlesInternalMP(&pix, fr, tr, pr_energy);
          break;
        case ptype::boundary:
          InteractParticlesBoundaryMP(&pix, fr, tr, pr_energy, virial);
          break;
        case ptype::tether:
          InteractParticlesTetherMP(&pix, fr, tr, pr_energy, virial);
          break;
        default:
          std::cout << "Wrong interaction type: " << (int)pix->type_ << std::endl;
          exit(1);
          break;
      }

      // Determine the interaction type, and run appropriate function
    } // omp for schedule(runtime) nowait (over interactions)
   
    // Reduce once all threads have finished
    #ifdef ENABLE_OPENMP
    #pragma omp barrier
    #endif
    ieh::ReduceMPRegion(tid, nsimples_, nspecies_, nthreads_, 
        &frc_, &trqc_, &virial_, &prc_energy_);

    delete[] fr;
    delete[] tr;
    delete[] virial;
  } // pragma omp parallel

  ReduceParticlesMP();
}

// External interactions
void InteractionEngine::InteractParticlesExternalMP(interaction_t **pix,
                                                      double **fr,
                                                      double **tr,
                                                      double *pe,
                                                      double **virial) {
  int idx = (*pix)->idx_; // actual locations
  int jdx = (*pix)->jdx_;
  auto part1 = (*simples_)[idx];
  auto part2 = (*simples_)[jdx];
  //if (part1->GetCID() == part2->GetCID()) return;

  // Calculate the potential, which we already have from the interaction
  PotentialBase *pot = (*pix)->pot_;
  if (pot->IsKMC()) return; // XXX FIXME is this necessary, guranteed to be the correct thing...

  // Minimum distance calc
  interactionmindist idm;
  MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);
  if ((idm.dr_mag2) > pot->GetRCut2()) return;

  // Already have manifest the location, no need to do an oid mapping
  auto sid2x = spec_ind_map_[part2->GetSID()];

  // Fire off the potential calc
  double fepot[4] = {0};
  pot->CalcPotential(&idm, part1, part2, fepot);

  #ifdef DEBUG
  if (debug_trace) {
    std::cout << "\tPOT EXTERNAL Interacting[" << idx << "," << part1->GetOID()
      << ":" << jdx << "," << part2->GetOID() << "] u: " << std::setprecision(16)
      << fepot[ndim_] << ", f: (" << fepot[0] << ", " << fepot[1];
    if (ndim_ == 3) {
      std::cout << ", " << fepot[2];
    }
    std::cout << ")\n";
  }
  #endif

  // Potential Energies
  pe[idx] += fepot[ndim_];
  pe[jdx] += fepot[ndim_];

  // Do the forces
  for (int i = 0; i < ndim_; ++i) {
      fr[i][idx] += fepot[i];
      fr[i][jdx] -= fepot[i];
      //Calculate virial only on particle two
      //FIXME This shouldn't be by species but by potential
      for(int j = i; j < ndim_; ++j)
        virial[3*i+j][sid2x] = virial[3*j+i][sid2x] += fr[i][jdx]*idm.dr[j];
  }

  // Calculate the torques
  double tau[3];
  cross_product(idm.contact1, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][idx] += tau[i];
  }
  cross_product(idm.contact2, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][jdx] -= tau[i];
  }
}

// KMC interactions
void InteractionEngine::InteractParticlesKMCMP(interaction_t **pix,
                                                 double **fr,
                                                 double **tr,
                                                 double *pe,
                                                 double **virial) {
  int idx = (*pix)->idx_;
  int jdx = (*pix)->jdx_;
  auto part1 = (*simples_)[idx];
  auto part2 = (*simples_)[jdx];

  auto neighbor = (*pix)->neighbor_;
  neighbor->kmc_= 0.0;
  // Check part1 for if we actually need to calculate or not, if not, just move on
  if (!part1->ApplyKMCInteraction()) {
    #ifdef DEBUG
    if (debug_trace) {
      std::cout << part1->GetOID() << " not applying KMC interaction\n";
    }
    #endif
    return;
  }

  // Calculate the potential
  PotentialBase *pot = (*pix)->pot_;
  #ifdef DEBUG
  if (!pot->IsKMC()) {
    std::cout << "InteractParticlesKMCMP wrong potential somehow, exiting\n";
    exit(1);
  }
  #endif

  // Check the neighbor list vs. interaction
  #ifdef DEBUG
  // The neighbor references the position in the LOCAL list of 
  // particles, not the global scope, so check to make sure that
  // they're the same thing
  if (neighbor->g_idx_ != (*pix)->jdx_) {
    std::cout << "Interaction neighbor translation wrong! Interaction: [" << (*pix)->idx_
      << ", " << (*pix)->jdx_ << "], global neighbor idx: " << neighbor->g_idx_
      << ", local neighbor idx: " << neighbor->idx_ << std::endl;
    exit(1);
  }
  #endif

  // Minimum distance calc
  interactionmindist idm;
  MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);
  if ((idm.dr_mag2) > pot->GetRCut2()) return;

  // Fire off the potential calc
  double fepot[4] = {0};
  pot->CalcPotential(&idm, part1, part2, fepot);

  #ifdef DEBUG
  if (debug_trace) {
    std::cout << "\tKMC Interacting[" << idx << "," << part1->GetOID()
      << ":" << jdx << "," << part2->GetOID() << "] kmc: " << std::setprecision(16)
      << fepot[ndim_] << ", (rmag2: " << idm.dr_mag2 << ")\n";
  }
  #endif

  neighbor->kmc_ = fepot[ndim_];
}

// Internal interactions
void InteractionEngine::InteractParticlesInternalMP(interaction_t **pix,
                                                      double **fr,
                                                      double **tr,
                                                      double *pe) {
  int idx = (*pix)->idx_;
  int jdx = (*pix)->jdx_;
  auto part1 = (*simples_)[idx];
  auto part2 = (*simples_)[jdx];

  // Check to see if we are doing this calculation (xlinks not being doubley bound...)
  bool do_calc = true;
  do_calc &= part1->ApplyInternalForce();
  do_calc &= part2->ApplyInternalForce();
  if (!do_calc) return;

  interactionmindist idm;
  MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);

  // Do the potential calc
  double fepot[4] = {0};
  PotentialBase *pot = (*pix)->pot_;
  pot->CalcPotential(&idm, part1, part2, fepot);

  #ifdef DEBUG
  if (debug_trace) {
    std::cout << "\tPOT INTERNAL Interacting[" << idx << "," << part1->GetOID()
      << ":" << jdx << "," << part2->GetOID() << "] u: " << std::setprecision(16)
      << fepot[ndim_] << ", f: (" << fepot[0] << ", " << fepot[1];
    if (ndim_ == 3) {
      std::cout << ", " << fepot[2];
    }
    std::cout << ")\n";
  }
  #endif

  // Potential Energies
  pe[idx] += fepot[ndim_];
  pe[jdx] += fepot[ndim_];

  // Do the forces
  for (int i = 0; i < ndim_; ++i) {
      fr[i][idx] += fepot[i];
      fr[i][jdx] -= fepot[i];
  }

  // Calculate the torques
  double tau[3];
  cross_product(idm.contact1, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][idx] += tau[i];
  }
  cross_product(idm.contact2, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][jdx] -= tau[i];
  }
}

// Boundary interactions
void InteractionEngine::InteractParticlesBoundaryMP(interaction_t **pix,
                                                      double **fr,
                                                      double **tr,
                                                      double *pe,
                                                      double **virial) {
  // Only 1 particle
  int idx = (*pix)->idx_;
  auto part1 = (*simples_)[idx];
  auto part2 = nullptr;

  // Always just do the calculation, boundary potentials know about their own
  // interaction distance calc
  interactionmindist idm;

  // Do the potential calc
  double fepot[4] = {0};
  PotentialBase *pot = (*pix)->pot_;
  pot->CalcPotential(&idm, part1, part2, fepot);

  #ifdef DEBUG
  if (debug_trace && fepot[ndim_] > 0.0) {
    std::cout << "\tPOT BOUNDARY Interacting[" << idx << "," << part1->GetOID()
      << "] u: " << std::setprecision(16)
      << fepot[ndim_] << ", f: (" << fepot[0] << ", " << fepot[1];
    if (ndim_ == 3) {
      std::cout << ", " << fepot[2];
    }
    std::cout << ")\n";
  }
  #endif

  // Potential Energies
  pe[idx] += fepot[ndim_];

  // Do the forces
  for (int i = 0; i < ndim_; ++i) {
      fr[i][idx] += fepot[i];
  }

  // Torques
  // Calculate the torques
  double tau[3];
  cross_product(idm.contact1, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][idx] += tau[i];
  }
}

// Tether interactions
void InteractionEngine::InteractParticlesTetherMP(interaction_t **pix,
                                                    double **fr,
                                                    double **tr,
                                                    double *pe,
                                                    double **virial) {
  // We cache the anchor information in the interaction in the case of
  // an anchor potential, so can just directly rip out of pix
  int idx = (*pix)->idx_;
  int jdx = (*pix)->jdx_;
  auto part1 = (*simples_)[idx];
  auto part2 = (*simples_)[jdx];
  auto sid2x = spec_ind_map_[part2->GetSID()];

  // Anchor information
  anchor_t *manchor = (*pix)->anchor_;

  // Calculating this potential is strange, since the minimum distance calculation is dependent
  // on the anchor point, and the potential tip, so call with the anchor point in the first position,
  // and part 2 in the second
  double rx0[3] = {0.0, 0.0, 0.0};
  double sx0[3] = {0.0, 0.0, 0.0};
  double rx1[3] = {0.0, 0.0, 0.0};
  double sx1[3] = {0.0, 0.0, 0.0};
  std::copy(manchor->pos0_, manchor->pos0_+3, rx0);
  std::copy(manchor->pos1_, manchor->pos1_+3, rx1);
  double dr[3] = {0.0, 0.0, 0.0};
  separation_vector(ndim_, nperiodic_, rx0, sx0, rx1, sx1, space_->unit_cell, dr);
  interactionmindist idm;
  std::copy(dr, dr+3, idm.dr);

  // Calculate the rcontacts for us, ugh (based on absolute vs relative positions)
  for (int i = 0; i < ndim_; ++i) {
    idm.contact1[i] = rx0[i] - part1->GetRigidPosition()[i];
    idm.contact2[i] = rx1[i] - part2->GetRigidPosition()[i];
  }

  // Fire off the potential calculation
  double fepot[4] = {0.0};
  PotentialBase *pot = (*pix)->pot_;
  pot->CalcPotential(&idm, part1, part2, fepot);

  #ifdef DEBUG
  if (debug_trace) {
    std::cout << "\tPOT TETHER Interacting[" << idx << "," << part1->GetOID()
      << ":" << jdx << "," << part2->GetOID() << "] u: " << std::setprecision(16)
      << fepot[ndim_] << ", f: (" << fepot[0] << ", " << fepot[1];
    if (ndim_ == 3) {
      std::cout << ", " << fepot[2];
    }
    std::cout << ")\n";
  }
  #endif

  // Potential Energies
  pe[idx] += fepot[ndim_];
  pe[jdx] += fepot[ndim_];

  // Do the forces
  for (int i = 0; i < ndim_; ++i) {
      fr[i][idx] += fepot[i];
      fr[i][jdx] -= fepot[i];
      //Calculate virial only on particle two
      //FIXME This shouldn't be by species but by potential
      for(int j = i; j < ndim_; ++j)
        virial[3*i+j][sid2x] = virial[3*j+i][sid2x] += fr[i][jdx]*idm.dr[j];
  }

  // Calculate the torques
  double tau[3];
  cross_product(idm.contact1, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][idx] += tau[i];
  }
  cross_product(idm.contact2, fepot, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
      tr[i][jdx] -= tau[i];
  }
}


// Reduce the particles and load the forces back onto them
void InteractionEngine::ReduceParticlesMP() {
  for (int i = 0; i < nsimples_; ++i) {
    auto part = (*simples_)[i];
    #ifdef DEBUG
    int oidx = (*oid_position_map_)[part->GetOID()];
    assert(i == oidx);
    #endif
    //if (i != oidx) {
    //  std::cout << "Isn't " << i << " = " << oidx << std::endl;
    //}
    double subforce[3] = {0.0, 0.0, 0.0};
    double subtorque[3] = {0.0, 0.0, 0.0};
    for (int idim = 0; idim < 3; ++idim) {
      subforce[idim] = frc_[idim*nsimples_+i];
      subtorque[idim] = trqc_[idim*nsimples_+i];
    }
    part->AddForceTorqueEnergy(subforce, subtorque, prc_energy_[i]);
  }

  //Add virial components to species
  for (int i = 0; i < nspecies_; ++i) {
    auto spec = (*species_)[i];
    double subvirial[9] = {};
    for (int j = 0; j < 9; ++j){
      subvirial[j] = virial_[j*nspecies_+i];
    }
    //std::copy(subvirial, subvirial+9, std::ostream_iterator<double>(std::cout,", "));
    spec->SetVirial(subvirial);
  }
}

// Dump functionality
void InteractionEngine::Dump() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "InteractionEngine::Dump\n";
  DumpSpecies();
  DumpInteractions();
  #endif
}

// Dump Species information
void InteractionEngine::DumpSpecies() {
  #ifdef DEBUG
  if (debug_trace) {
    std::cout << "----------------\n";
    std::cout << "DumpSpecies\n";
    for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
      std::cout << "Species[" << (int)(*spec)->GetSID() << "] ->\n";
      (*spec)->Dump();
    }
  }
  #endif
}

// Dump interaction information
void InteractionEngine::DumpInteractions() {
  #ifdef DEBUG
  std::cout << "----------------\n";
  std::cout << "DumpInteractions\n";
  for (int ixs = 0; ixs < interactions_->size(); ++ixs) {
    auto mixs = (*interactions_)[ixs];
    std::cout << "[" << ixs << "] {" << mixs.idx_ << " -> " << mixs.jdx_ << "}, type: "
      << PtypeToString(mixs.type_);
    if (mixs.kmc_target_ != SID::none) {
      std::cout << ", kmc_target: " << SIDToString(mixs.kmc_target_)
        << std::setprecision(16) << ", kmc: " << mixs.neighbor_->kmc_;
    }
    std::cout << std::endl;
  }
  #endif
}
