#include "br_rod.h"

void BrRod::Init() {
  Simple::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 4*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  orientation_[0] = 0; // Set default color to red
}
void BrRod::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  Integrate();
  UpdatePeriodic();
  ClearInteractions();
}

// Basic verlet integrator, very stable
void BrRod::Integrate() {
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

void BrRodSpecies::InitPotentials (system_parameters *params) {
  AddPotential(SID::br_rod, SID::br_rod, 
      // Set br_rod-br_rod interaction
      new LJ126(params->lj_epsilon,params->rod_diameter,
        space_, pow(2, 1.0/6.0)*params->br_rod_diameter));
}
