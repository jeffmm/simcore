#include "simple.h"

unsigned int Simple::next_oid_ = 0;

Simple::Simple() {
  oid_ = ++next_oid_;
}

Simple::Simple(int n_dim) {
  oid_ = ++next_oid_;
  Init(n_dim);
}

Simple::~Simple() {
  delete position_;
  delete orientation_;
  delete velocity_;
  delete force_;
}

void Simple::Init(int n_dim) {
  n_dim_ = n_dim;
  position_ = new double(n_dim_);
  orientation_ = new double(n_dim_);
  velocity_ = new double(n_dim_);
  force_ = new double(n_dim_);
  radius_ = 1;
  length_ = 0;
  ZeroForce();
}

const unsigned int Simple::GetSimpleID() const {
  return oid_;
}

void Simple::SetRadius(double new_radius) {
  radius_ = new_radius;
}

void Simple::SetLength(double new_length) {
  length_ = new_length;
}

double Simple::GetRadius() const {
  return radius_;
}

double Simple::GetLength() const {
  return length_;
}

void Simple::ZeroForce() {
  for (int i=0; i<n_dim_; ++i)
    force_[i] = 0;
}

void Simple::SetPosition(double *new_position) {
  for (int i=0; i<n_dim_; ++i)
    position_[i] = new_position[i];
}

void Simple::SetOrientation(double *new_orientation) {
  for (int i=0; i<n_dim_; ++i)
    orientation_[i] = new_orientation[i];
}

void Simple::SetVelocity(double *new_velocity) {
  for (int i=0; i<n_dim_; ++i)
    velocity_[i] = new_velocity[i];
}

double *Simple::GetPosition() {
  return position_;
}

double *Simple::GetOrientation() {
  return orientation_;
}

double *Simple::GetVelocity() {
  return velocity_;
}

double *Simple::GetForce() {
  return force_;
}

void Simple::AddInteraction(interaction new_interaction) {
  interactions.push_back(new_interaction);
}

void Simple::ApplyInteractions() {}

