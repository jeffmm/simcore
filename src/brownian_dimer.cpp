#include "brownian_dimer.h"

void BrownianDimer::InitRandom(double sys_radius) {
  generate_random_unit_vector(n_dim_, com_, rng_.r);
  for (int i=0; i<n_dim_; ++i) 
    com_[i] *= gsl_rng_uniform_pos(rng_.r)*sys_radius;
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  double z = 0.5*length_;
  double pos[3];
  for (auto i_bead = beads_.begin(); i_bead != beads_.end(); ++i_bead) {
    for (int i=0; i<n_dim_; ++i)
      pos[i] = com_[i] + z*orientation_[i];
    i_bead->SetPosition(pos);
    i_bead->SetDiameter(diameter_);
    z-=length_;
  }
  SetDiffusion();
  UpdateOrientation();
}

void BrownianDimer::InitRest() {
  double z = 0.5*length_;
  double *pos = new double[n_dim_];
  for (int i=0; i<n_dim_; ++i)
    pos[i] = 0;
  for (auto i_bead = beads_.begin(); i_bead != beads_.end(); ++i_bead) {
    pos[n_dim_-1] = z;
    i_bead->SetPosition(pos);
    i_bead->SetDiameter(diameter_);
    z-=length_;
  }
  delete[] pos;
  SetDiffusion();
  UpdateOrientation();
}

void BrownianDimer::KickBeads() {
  for (int i_bead=0; i_bead<2; ++i_bead) {
    for (int i=0; i<n_dim_; ++i) {
      double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
      force_[i] = kick*diffusion_;
    }
    beads_[i_bead].AddForce(force_);
  }
}

void BrownianDimer::UpdateOrientation() {
  double const * const r1=beads_[0].GetPosition();
  double const * const r2=beads_[1].GetPosition();
  length_ = 0;
  for (int i=0; i<n_dim_; ++i)
    length_ += SQR(r2[i]-r1[i]);
  length_=sqrt(length_);
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = (r2[i]-r1[i])/length_;
    com_[i] = r1[i] + 0.5*length_*orientation_[i];
  }
  beads_[0].SetOrientation(orientation_);
  beads_[1].SetOrientation(orientation_);
}

void BrownianDimer::InternalForces() {
  for (int i=0; i<n_dim_; ++i)
    force_[i] = k_spring_ * (length_-eq_length_) * orientation_[i];
  beads_[0].AddForce(force_);
  for (int i=0; i<n_dim_; ++i)
    force_[i] = -force_[i];
  beads_[1].AddForce(force_);
}

void BrownianDimer::UpdatePosition() {
  ZeroForces();
  KickBeads();
  InternalForces();
  Integrate();
}

void BrownianDimer::ZeroForces() {
  for (int i_bead=0; i_bead<2; ++i_bead)
    beads_[i_bead].ZeroForce();
}

void BrownianDimer::Integrate() {
  double pos[3] = {0, 0, 0};
  for (int i_bead=0; i_bead<2; ++i_bead) {
    double const * const r = beads_[i_bead].GetPosition();
    double const * const f = beads_[i_bead].GetForce();
    for (int i=0; i<n_dim_; ++i)
      pos[i] = r[i] + f[i] * delta_ / diameter_;
    beads_[i_bead].SetPosition(pos);
  }
  UpdateOrientation();
}

void BrownianDimer::Draw(std::vector<graph_struct*> * graph_array) {
  UpdateOrientation();
  double const * const r1 = beads_[0].GetPosition();
  memcpy(g_.r,com_,sizeof(com_));
  memcpy(g_.u,orientation_,sizeof(orientation_));
  g_.length = length_;
  g_.diameter = 0.5*diameter_;
  graph_array->push_back(&g_);
  for (int i=0; i<2; ++i)
    beads_[i].Draw(graph_array);
}


