#include "br_bead.hpp"

BrBead::BrBead() : Object() {}

void BrBead::Init(br_bead_parameters *sparams) {
  sparams_ = sparams;
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  driving_factor_ = sparams_->driving_factor;
  stoch_flag_ = params_->stoch_flag; // flag for thermal forces
  SetDiffusion();
  InsertBrBead();
}

int BrBead::GetCount() { return 1; }

void BrBead::ZeroForce() { Object::ZeroForce(); }

void BrBead::InsertBrBead() {
  if (sparams_->insertion_type.compare("random") == 0) {
    InsertRandom();
  } else if (sparams_->insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (sparams_->insertion_type.compare("centered_random") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  } else if (sparams_->insertion_type.compare("centered_oriented") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (sparams_->insertion_type.compare("custom") == 0) {
    // Nothing to do
  } else {
    Logger::Error("BrBead insertion type not recognized!");
  }
}

void BrBead::UpdatePosition() {
  SetPrevPosition(position_);
  ApplyForcesTorques();
  Integrate();
  UpdatePeriodic();
}

void BrBead::ApplyForcesTorques() {
  // Add random thermal kick to the bead
  if (stoch_flag_) {
    for (int i = 0; i < n_dim_; ++i) {
      double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
      force_[i] += kick * diffusion_;
    }
  }
  if (driving_factor_ > 0) {
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] += driving_factor_ * orientation_[i];
    }
  }
}

void BrBead::SetDiffusion() {
  gamma_trans_ = 1.0 / (diameter_);
  gamma_rot_ = 3.0 / CUBE(diameter_);
  diffusion_ = sqrt(24.0 * diameter_ / delta_);
}

void BrBead::Translate() {
  double dr[3];
  double fmag = 0;
  for (int i = 0; i < n_dim_; ++i) {
    fmag += force_[i] * force_[i];
    dr[i] = force_[i] * delta_ * gamma_trans_;
    position_[i] += dr[i];
  }
}

void BrBead::Rotate() {
  double unit_torque[3], temp[3], r_rel[3];
  double domega, cos_domega, sin_domega, torque_mag;
  // First rotate orientation vector of sphere
  if (n_dim_ == 2) {
    domega = torque_[2] * delta_ * gamma_rot_;
    cos_domega = cos(domega);
    sin_domega = sin(domega);
    std::copy(orientation_, orientation_ + 3, temp);
    orientation_[0] = cos_domega * temp[0] - sin_domega * temp[1];
    orientation_[1] = sin_domega * temp[0] + cos_domega * temp[1];
  } else if (n_dim_ == 3) {
    torque_mag = 0.0;
    for (int i = 0; i < 3; ++i)
      torque_mag += torque_[i];
    for (int i = 0; i < 3; ++i)
      unit_torque[i] = torque_[i] / torque_mag;
    domega = torque_mag * delta_ * gamma_rot_;
    rotate_vector(orientation_, unit_torque, domega);
  }
  normalize_vector(orientation_, n_dim_);
}

void BrBead::Integrate() {
  Translate();
  // Rotate();
}

void BrBead::GetInteractors(std::vector<Object *> *ix) { ix->push_back(this); }

void BrBead::Draw(std::vector<graph_struct *> *graph_array) {
  Object::Draw(graph_array);
}
