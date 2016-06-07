#include "md_bead.h"

void MDBead::Init() {
  Bead::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 4*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  orientation_[0] = 0; // Set default color to red
}
void MDBead::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  Integrate();
  UpdatePeriodic();
  ClearInteractions();
}

// Basic verlet integrator, very stable
void MDBead::Integrate() {
  double delta2 = SQR(delta_);
  double temp_pos[3];
  for (int i=0; i<n_dim_; ++i) {
    temp_pos[i] = position_[i];
    position_[i] = 2.0*position_[i] - prev_position_[i] 
      + ((force_[i]/mass_))*delta2;
    velocity_[i] =  (position_[i] - prev_position_[i])/(2.0*delta_);
    prev_position_[i] = temp_pos[i];
  }
}

void MDBead::UpdateKineticEnergy() {
  double vel_mag_sqr = 0.0;
  for (int i=0; i<n_dim_; ++i)
    vel_mag_sqr += SQR(velocity_[i]);
  k_energy_ = 0.5 * mass_ * vel_mag_sqr;
}

double const MDBead::GetKineticEnergy() {
  UpdateKineticEnergy();
  return k_energy_;
}

void MDBeadSpecies::InitPotentials (system_parameters *params) {
  AddPotential(SID::md_bead, SID::md_bead, 
      // Set md_bead-md_bead interaction
      new LJ126(params->lj_epsilon,params->md_bead_diameter,
        space_, 2.5*params->md_bead_diameter));
}
