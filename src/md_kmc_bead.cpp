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
    Integrate();
    UpdatePeriodic();
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
void MDKMCBead::UpdateProbability() {
  n_exp_ = 0.0;
  double binding_affinity = eps_eff_ * on_rate_ * delta_;
  n_exp_ += binding_affinity * kmc_energy_;
}

void MDKMCBead::StepKMC() {
  if (!bound_ && n_exp_ > 0.0) {
    double roll = gsl_rng_uniform(rng_.r);
    if (roll < n_exp_) {
      printf("Successful KMC move {nexp: %2.2f}, {roll: %2.2f}\n", n_exp_, roll);
      bound_ = true;
      // Figure out which friend to fall on, whooooos it gonna be....
    }
  }
}


// Species kmc stuff
void MDKMCBeadSpecies::PrepKMC() {
  UpdateProbability();
}

void MDKMCBeadSpecies::UpdateProbability() {
  n_exp_ = 0.0;
  for (auto kmcit = members_.begin(); kmcit != members_.end(); ++kmcit) {
    (*kmcit)->UpdateProbability();
    n_exp_ += (*kmcit)->GetNExp();
  }
}

void MDKMCBeadSpecies::StepKMC() {
  for (auto kmcit = members_.begin(); kmcit != members_.end(); ++kmcit) {
    (*kmcit)->StepKMC();
  }
}
