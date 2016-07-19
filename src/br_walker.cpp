#include "br_walker.h"

void BrWalker::Init() {
  Simple::Init();
  orientation_[0] = 0; // Set default color to red
  orientation_[1] -=orientation_[1];
  n_exp_ = 0.0;
}

void BrWalker::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}

void BrWalker::UpdatePosition() {
  KickBead();
  ApplyInteractions();
  for (int i=0; i<n_dim_; ++i)
    position_[i] = position_[i] + force_[i] * delta_ / diameter_;
  UpdatePeriodic();
  ClearInteractions();
  ZeroForce();
}

void BrWalker::UpdatePositionMP() {
  // If bound, the position is set by the kmc engine module
  if (!bound_) {
    KickBead();
    for (int i = 0; i < n_dim_; ++i) {
        position_[i] = position_[i] + force_[i] * delta_ / diameter_;
        dr_tot_[i] += position_[i] + force_[i] * delta_ / diameter_;
    }
  } else {
    for (int i = 0; i < n_dim_; ++i) {
      dr_tot_[i] += position_[i] - prev_position_[i];
    }
  }
  UpdatePeriodic();
}

// kmc stuff
void BrWalker::PrepKMC(std::vector<neighbor_t>* neighbors) {

}

void BrWalker::StepKMC() {

}

void BrWalker::Attach(int idx, double pos) {
  attachidx_ = idx;
  attachpos_ = pos;
}
