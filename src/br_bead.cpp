#include "simcore/br_bead.hpp"

BrBead::BrBead(unsigned long seed) : Object(seed) {
  SetSID(species_id::br_bead);
}

void BrBead::Init(br_bead_parameters *sparams) {
  sparams_ = sparams;
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  driving_factor_ = sparams_->driving_factor;
  zero_temperature_ = params_->zero_temperature;
  chiral_handedness_ = sparams_->chiral_handedness;
  alignment_interaction_ = sparams_->alignment_interaction;
  alignment_torque_ = sparams_->alignment_torque;
  noise_factor_ = sparams_->noise_factor;
  if (sparams_->randomize_handedness) {
    chiral_handedness_ = (rng_.RandomUniform() > 0.5 ? 1 : -1);
  }
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (chiral_handedness_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
  driving_torque_ = sparams_->driving_torque;
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
    rng_.RandomUnitVector(n_dim_, orientation_);
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
  SetPrevPosition(GetPosition());
  SetPrevOrientation(GetOrientation());
  ApplyForcesTorques();
  Integrate();
  UpdatePeriodic();
}

void BrBead::ApplyForcesTorques() {
  // Add random thermal kick to the bead
  if (!zero_temperature_) {
    for (int i = 0; i < n_dim_; ++i) {
      double kick = rng_.RandomUniform() - 0.5;
      force_[i] += kick * diffusion_;
    }
    if (n_dim_ == 2) {
      double kick = rng_.RandomUniform() - 0.5;
      torque_[2] += diffusion_rot_ * kick;
    }
  }
  if (driving_factor_ > 0) {
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] += driving_factor_ * orientation_[i];
    }
  }
  if (driving_torque_ > 0) {
    if (n_dim_ == 2) {
      torque_[2] += chiral_handedness_ * driving_torque_;
    }
  }
  if (alignment_interaction_) {
    for (auto ix = ixs_.begin(); ix != ixs_.end(); ++ix) {
      if (ix->first->pause_interaction)
        continue;
      double dr2 = ix->first->dr_mag2;
      double cp[3] = {0, 0, 0};
      double align_t = alignment_torque_;
      if (ix->second) {
        const double * const u = ix->first->obj2->GetOrientation();
        cross_product(orientation_, u, cp, n_dim_);
        int sign = SIGNOF(cp[2]);
        align_t *= sign*0.5*(1 - dot_product(n_dim_, u, orientation_));
      } else {
        const double * const u = ix->first->obj1->GetOrientation();
        cross_product(orientation_, u, cp, n_dim_);
        int sign = SIGNOF(cp[2]);
        align_t *= sign*0.5*(1 - dot_product(n_dim_, u, orientation_));
      }
      cp[2] = align_t/dr2;
      AddTorque(cp);
    }
  }
}

void BrBead::SetDiffusion() {
  gamma_trans_ = 1.0 / (diameter_);
  gamma_rot_ = 3.0 / CUBE(diameter_);
  diffusion_ = noise_factor_ * sqrt(24.0 * diameter_ / delta_);
  diffusion_rot_ = noise_factor_ * sqrt(8.0 * CUBE(diameter_) / delta_);
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
    rotate_vector(orientation_, unit_torque, domega, n_dim_);
  }
  normalize_vector(orientation_, n_dim_);
}

void BrBead::Integrate() {
  Translate();
  Rotate();
}

void BrBead::GetInteractors(std::vector<Object *> &ix) { ix.push_back(this); }

void BrBead::Draw(std::vector<graph_struct *> &graph_array) {
  Object::Draw(graph_array);
}
void BrBead::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing br_bead specs, object id: %d", GetOID());
  Object::WriteSpec(ospec);
  ospec.write(reinterpret_cast<char *>(&chiral_handedness_), sizeof(int));
}

void BrBead::ReadSpec(std::fstream &ispec) {
  Object::ReadSpec(ispec);
  ispec.read(reinterpret_cast<char *>(&chiral_handedness_), sizeof(int));
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (chiral_handedness_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
}


