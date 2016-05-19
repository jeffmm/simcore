#include "brownian_dimer.h"

void BrownianDimer::InitRandom(double sys_radius) {
  generate_random_unit_vector(n_dim_, position_, rng_.r);
  for (int i=0; i<n_dim_; ++i) 
    position_[i] *= gsl_rng_uniform_pos(rng_.r)*sys_radius;
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  double z = 0.5*length_;
  double pos[3];
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    for (int i=0; i<n_dim_; ++i)
      pos[i] = position_[i] + z*orientation_[i];
    i_bead->SetPosition(pos);
    i_bead->SetDiameter(diameter_);
    z-=length_;
  }
  UpdatePeriodic();
  UpdateOrientation();
}

void BrownianDimer::KickBeads() {
  for (auto i_bead=elements_.begin(); i_bead!= elements_.end(); ++i_bead)
    i_bead->KickBead();
}

void BrownianDimer::UpdateOrientation() {
  double const * const r1=elements_[0].GetPosition();
  double const * const r2=elements_[1].GetPosition();
  double const * const s1=elements_[0].GetScaledPosition();
  double const * const s2=elements_[1].GetScaledPosition();
  length_ = 0;
  double dr[3];
  separation_vector(n_dim_, space_->n_periodic, r1, s1, r2, s2, space_->unit_cell, dr);
  for (int i=0; i<n_dim_; ++i) {
    length_ += SQR(dr[i]);
  }
  length_=sqrt(length_);
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = (dr[i])/length_;
    position_[i] = r1[i] + 0.5*length_*orientation_[i];
  }
  for (auto i_bead=elements_.begin(); i_bead!= elements_.end(); ++i_bead)
    i_bead->SetOrientation(orientation_);
}

void BrownianDimer::InternalForces() {
  for (int i=0; i<n_dim_; ++i)
    force_[i] = k_spring_ * (length_-eq_length_) * orientation_[i];
  elements_[0].AddForce(force_);
  for (int i=0; i<n_dim_; ++i)
    force_[i] = -force_[i];
  elements_[1].AddForce(force_);
}

void BrownianDimer::UpdatePosition() {
  ZeroForces();
  KickBeads();
  InternalForces();
  Integrate();
}

void BrownianDimer::Integrate() {
  double pos[3] = {0, 0, 0};
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    double const * const r = i_bead->GetPosition();
    double const * const f = i_bead->GetForce();
    for (int i=0; i<n_dim_; ++i)
      pos[i] = r[i] + f[i] * delta_ / diameter_;
    i_bead->SetPosition(pos);
  }
  UpdatePeriodic();
  UpdateOrientation();
}

void BrownianDimer::Draw(std::vector<graph_struct*> * graph_array) {
  // draw half of rod coming from each bead (looks good for periodic boundaries)
  double const * const r1 = elements_[0].GetPosition();
  double const * const r2 = elements_[1].GetPosition();
  double r[3];
  for (int i=0; i<n_dim_; ++i)
    r[i] = r1[i] + 0.25*length_*orientation_[i];
  memcpy(g_.r,r,sizeof(r));
  memcpy(g_.u,orientation_,sizeof(orientation_));
  g_.length = 0.5 * length_;
  g_.diameter = 0.5*diameter_;
  graph_array->push_back(&g_);
  for (int i=0; i<n_dim_; ++i)
    r[i] = r2[i] - 0.25*length_*orientation_[i];
  memcpy(g2_.r,r,sizeof(r));
  memcpy(g2_.u,orientation_,sizeof(orientation_));
  g2_.length = 0.5 * length_;
  g2_.diameter = 0.5*diameter_;
  graph_array->push_back(&g2_);
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead)
    i_bead->Draw(graph_array);
}


