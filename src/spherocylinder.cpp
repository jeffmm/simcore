#include "simcore/spherocylinder.hpp"

Spherocylinder::Spherocylinder(unsigned long seed) : BrRod(seed) {
  SetSID(species_id::spherocylinder);
}

void Spherocylinder::SetParameters() {
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  length_ = sparams_->length;
  std::fill(body_frame_, body_frame_ + 6, 0.0);
  SetDiffusion();
}

void Spherocylinder::Init(spherocylinder_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  InsertSpherocylinder();
  interactors_.push_back(this);
}

void Spherocylinder::InsertSpherocylinder() {
  InsertRod(sparams_->insertion_type);
}

void Spherocylinder::UpdatePosition() {
  SetPrevPosition(position_);
  ApplyForcesTorques();
  Integrate();
  UpdatePeriodic();
}

// TODO: Interacting spheros
void Spherocylinder::ApplyForcesTorques() {}

