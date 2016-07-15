#include "md_kmc_bead.h"

void MDKMCBead::Init() {
  Simple::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 4*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  orientation_[0] = 0; // Set default color to red
  orientation_[1] -=orientation_[1];
  n_exp_ = 0.0;
}
void MDKMCBead::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  Integrate();
  UpdatePeriodic();
  ClearInteractions();
}

// Update the position based on MP
void MDKMCBead::UpdatePositionMP() {
  if (!bound_) {
    Integrate();
  } else {
    SelfSetPosition();
  }
  UpdatePeriodic();
}

// We are attached to someone else somewhere, so use
// that instead, just calculate th dr_tot
void MDKMCBead::SelfSetPosition() {
  for (int i = 0; i < n_dim_; ++i) {
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
}

// Basic verlet integrator, very stable
void MDKMCBead::Integrate() {
  double delta2 = SQR(delta_);
  double temp_pos[3];
  for (int i=0; i<n_dim_; ++i) {
    temp_pos[i] = position_[i];
    position_[i] = 2.0*position_[i] - prev_position_[i] 
      + ((force_[i]/mass_))*delta2;
    velocity_[i] =  (position_[i] - prev_position_[i])/(2.0*delta_);
    prev_position_[i] = temp_pos[i];
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
}

void MDKMCBead::UpdateKineticEnergy() {
  double vel_mag_sqr = 0.0;
  for (int i=0; i<n_dim_; ++i)
    vel_mag_sqr += SQR(velocity_[i]);
  k_energy_ = 0.5 * mass_ * vel_mag_sqr;
}

double const MDKMCBead::GetKineticEnergy() {
  UpdateKineticEnergy();
  return k_energy_;
}

// kmc stuff
void MDKMCBead::PrepKMC(std::vector<neighbor_t>* neighbors) {
  n_exp_ = 0.0;
  if (bound_) return; //NOT expected to bind in this timestep
  neighbors_ = neighbors;
  double binding_affinity = eps_eff_ * on_rate_ * delta_;
  for (auto nldx = neighbors_->begin(); nldx != neighbors_->end(); ++nldx) {
    n_exp_ += binding_affinity * nldx->kmc_;
  }
}

void MDKMCBead::StepKMC() {
}

void MDKMCBead::Attach(int idx) {
  attachidx_ = idx;
}


// Species kmc stuff
void MDKMCBeadSpecies::PrepKMC() {
  n_exp_ = 0.0;
  for (auto kmcit = members_.begin(); kmcit != members_.end(); ++kmcit) {
    n_exp_ += (*kmcit)->GetNExp();
  }
}

void MDKMCBeadSpecies::StepKMC() {
  nbound_ = 0;
  nfree_ = 0;
  for (auto kmcit = members_.begin(); kmcit != members_.end(); ++kmcit) {
    if ((*kmcit)->GetBound()) {
      nbound_++;
    } else {
      nfree_++;
    }
    (*kmcit)->StepKMC();
  }
}

void MDKMCBeadSpecies::DumpKMC() {
  // print out the information appropriate to kmc
  if (debug_trace) {
    printf("Species[MDKMCBeadSpecies] -> dump\n");
    printf("\t{n_exp(this delta): %2.4f}\n", n_exp_);
    printf("\t{nfree:  %d}\n", nfree_);
    printf("\t{nbound: %d}\n", nbound_);
    for (auto kmcit = members_.begin(); kmcit != members_.end(); ++kmcit) {
      printf("\t\t[%d] -> {n_exp: %2.4f}, {bound: %s}, {parent index: %d}\n", (*kmcit)->GetOID(), (*kmcit)->GetNExp(), (*kmcit)->GetBound() ? "true " : "false",
          (*kmcit)->GetAttach());
    }
  }
}
