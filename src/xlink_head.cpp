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
    std::copy(position_, position_+n_dim_, prev_position_);
    for (int i = 0; i < n_dim_; ++i) {
        position_[i] = position_[i] + force_[i] * delta_ / diameter_;
    }
  }
  for (int i = 0; i < n_dim_; ++i) {
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
  UpdatePeriodic();
}

void XlinkHead::Draw(std::vector<graph_struct*> * graph_array) {
  std::copy(position_, position_+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  std::copy(color_, color_+4, g_.color);
  g_.length = length_;
  g_.diameter = diameter_ + 0.5;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
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

void XlinkHead::Bind(int idx, int ridx, int cidx, double pos) {
  attachidx_ = idx;
  attachridx_ = ridx;
  attachcidx_ = cidx;
  attachpos_ = pos;
  bound_ = true;
}

bool XlinkHead::GetBind(int *idx, int *ridx, int *cidx, double *pos) {
  (*idx) = attachidx_;
  (*ridx) = attachridx_;
  (*cidx) = attachcidx_;
  (*pos) = attachpos_;
  return bound_;
}
