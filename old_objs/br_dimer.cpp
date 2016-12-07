#include "br_dimer.h"

void BrDimer::Init() {
  Composite::InsertRandom(length_+diameter_);
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  double z = 0.5*length_;
  double pos[3];
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    for (int i=0; i<n_dim_; ++i)
      pos[i] = position_[i] + z*orientation_[i];
    i_bead->SetPosition(pos);
    i_bead->SetDiameter(diameter_);
    i_bead->SetDiffusion();
    i_bead->SetSpace(space_);
    i_bead->UpdatePeriodic();
    z-=length_;
  }
  UpdateOrientation();
}

void BrDimer::UpdateOrientation() {
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
  UpdatePeriodic();
  for (auto i_bead=elements_.begin(); i_bead!= elements_.end(); ++i_bead)
    i_bead->SetOrientation(orientation_);
}

void BrDimer::InternalForces() {
  for (int i=0; i<n_dim_; ++i)
    force_[i] = k_spring_ * (length_-eq_length_) * orientation_[i];
  elements_[0].AddForce(force_);
  for (int i=0; i<n_dim_; ++i)
    force_[i] = -force_[i];
  elements_[1].AddForce(force_);
}

void BrDimer::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  InternalForces();
  Integrate();
}

void BrDimer::UpdatePositionMP() {
  InternalForces();
  Integrate();
}

void BrDimer::ApplyInteractions() {
   //We have no internal constraints
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->ApplyInteractions();
  }
}

void BrDimer::Integrate() {
  double pos[3] = {0, 0, 0};
  // the beads know how to update position and periodicity
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->UpdatePositionMP();
  }
  UpdateOrientation();
}

void BrDimer::Draw(std::vector<graph_struct*> * graph_array) {
  // draw half of rod coming from each bead (looks good for periodic boundaries)
  double const * const r1 = elements_[0].GetPosition();
  double const * const r2 = elements_[1].GetPosition();
  double r[3];
  for (int i=0; i<n_dim_; ++i)
    r[i] = r1[i] + 0.25*length_*orientation_[i];
  std::copy(r, r+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  g_.length = 0.5 * length_;
  g_.diameter = 0.2*diameter_;
  graph_array->push_back(&g_);
  for (int i=0; i<n_dim_; ++i)
    r[i] = r2[i] - 0.25*length_*orientation_[i];
  std::copy(r, r+3, g2_.r);
  std::copy(orientation_, orientation_+3, g2_.u);
  g2_.length = 0.5 * length_;
  g2_.diameter = 0.2*diameter_;
  graph_array->push_back(&g2_);
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead)
    i_bead->Draw(graph_array);
}

//void BrDimerSpecies::InitPotentials(system_parameters *params) {
  //AddPotential(SID::br_dimer, SID::br_dimer, new LJ126(params->lj_epsilon, params->br_bead_diameter, space_, 2.5*params->br_bead_diameter));
  //double dimer_neon_diam = 0.5*(params->br_bead_diameter + params->neon_diameter);
  //double dimer_argon_diam = 0.5*(params->br_bead_diameter + params->argon_diameter);
  //AddPotential(SID::neon, SID::br_dimer, new LJ126(params->lj_epsilon, dimer_neon_diam, space_, 2.5*dimer_neon_diam));
  //AddPotential(SID::argon, SID::br_dimer, new LJ126(params->lj_epsilon, dimer_argon_diam, space_, 2.5*dimer_argon_diam));
//}

