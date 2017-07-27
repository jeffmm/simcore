#include "bond.h"

void Bond::Init(Site *s1, Site *s2) {
  double const * const r1 = s1->GetPosition();
  double const * const r2 = s2->GetPosition();
  diameter_ = s1->GetDiameter();
  length_ = 0;
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = r2[i] - r1[i];
    length_ += orientation_[i]*orientation_[i];
  }
  length_ = sqrt(length_);
  for (int i=0;i<n_dim_; ++i) {
    orientation_[i] /= length_;
    position_[i] = r1[i] + 0.5* length_ * orientation_[i];
  }
  UpdatePeriodic();
  SetRigidPosition(GetPosition());
  SetRigidScaledPosition(GetScaledPosition());
  SetRigidOrientation(orientation_);
  SetRigidLength(length_);
  SetRigidDiameter(diameter_);
  // Set neighbor IDs
  SetNIDS(s1->GetOID(), s2->GetOID());
  SetNeighbors(s1, s2);
}

void Bond::Draw(std::vector<graph_struct*> * graph_array) {
  for (int i=0; i<space_->n_periodic; ++i)
    g_.r[i] = scaled_position_[i];
  for (int i=space_->n_periodic; i<n_dim_; ++i)
    g_.r[i] = position_[i];
  std::copy(orientation_, orientation_+3, g_.u);
  g_.color = color_;
  //g_.length = (length_-diameter_ > 0 ? length_-diameter_ : 0);
  g_.length = length_;
  if (graph_diameter_ == 0) 
    g_.diameter = diameter_;
  else 
    g_.diameter = graph_diameter_;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
}

