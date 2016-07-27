#include "xlink_head.h"

void XlinkHead::Init() {
  Simple::Init();
  orientation_[0] = 0; // Set default color to red
  orientation_[1] = 1;
  orientation_[2] = 1;
  n_exp_0_1_ = 0.0;
  n_exp_1_2_ = 0.0;
}

void XlinkHead::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}

void XlinkHead::UpdatePositionMP() {
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
void XlinkHead::PrepKMC(std::vector<neighbor_t>* neighbors) {

}

void XlinkHead::StepKMC() {

}

void XlinkHead::Attach(int idx, double pos) {
  attachidx_ = idx;
  attachpos_ = pos;
}
