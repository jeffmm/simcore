#include "md_bead.h"

void MDBead::Init() {
  Bead::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 4*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_force_[i] = force_[i];
  }
}
void MDBead::UpdatePosition() {
  Integrate();
  UpdatePeriodic();
  ZeroForce();
}
void MDBead::Integrate() {
  for (int i=0; i<n_dim_; ++i) {
    position_[i] = position_[i] + delta_ * velocity_[i] + 0.5* prev_force_[i] * SQR(delta_) / mass_;
    velocity_[i] = velocity_[i] + 0.5 * (prev_force_[i] + force_[i]) * delta_ / mass_;
    prev_force_[i] = force_[i];
  }
}
void MDBead::UpdateEnergy() {
  double vel_mag_sqr = 0.0;
  for (int i=0; i<n_dim_; ++i)
    vel_mag_sqr += SQR(velocity_[i]);
  energy_ = 0.5 * mass_ * vel_mag_sqr;
}
double const MDBead::GetEnergy() {
  UpdateEnergy();
  return energy_;
}

