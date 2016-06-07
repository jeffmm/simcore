#include "brownian_bead.h"

void BrownianBead::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}
void BrownianBead::UpdatePosition() {
  KickBead();
  ApplyInteractions();
  for (int i=0; i<n_dim_; ++i)
    position_[i] = position_[i] + force_[i] * delta_ / diameter_;
  UpdatePeriodic();
  ClearInteractions();
  ZeroForce();
}
void BrownianBeadSpecies::InitPotentials (system_parameters *params) {
  AddPotential(SID::brownian_bead, SID::brownian_bead, 
      new LJ126(params->lj_epsilon,params->br_bead_diameter,
                space_, 2.5*params->br_bead_diameter));
}


