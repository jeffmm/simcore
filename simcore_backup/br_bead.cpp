#include "br_bead.h"

void BrBead::Init() {
  Simple::Init();
  std::fill(body_frame_,body_frame_+6,0);
}

void BrBead::UpdatePosition() {
  Translate();
  Rotate();
  UpdatePeriodic();
}

void BrBead::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}

void BrBead::Translate() {
  // Add random component to force
  KickBead();
  Driving();
  std::copy(position_, position_+n_dim_, prev_position_);
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = position_[i] + force_[i] * delta_ / diameter_;
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
}

void BrBead::Driving() {
  double f_dr[3];
  for (int i=0; i<n_dim_; ++i)
    f_dr[i] = orientation_[i]*driving_factor_;
  AddForce(f_dr);
}

void BrBead::SetDiffusion() {
  friction_rot_ = 3.0*diameter_*diameter_*diameter_;
  rand_sigma_rot_ = sqrt(2.0*delta_/friction_rot_);
  diffusion_ = sqrt(24.0*diameter_/delta_);
}

void BrBead::Rotate() {
// Now handle the random orientation update
  GetBodyFrame();
  for (int j=0; j<n_dim_-1; ++j) {
    double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
    for (int i=0; i<n_dim_; ++i)
      orientation_[i] += mag * body_frame_[n_dim_*j+i];
  }
  normalize_vector(orientation_, n_dim_);
}

/* calculates vector(s) orthogonal to orientation */
void BrBead::GetBodyFrame() {
  if (n_dim_==2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  }
  else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]),n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

