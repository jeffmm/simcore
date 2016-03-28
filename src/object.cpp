#include "object.h"

unsigned int Object::next_oid_ = 0;

Object::Object() {
  oid_ = ++next_oid_;
}

Object::Object(int n_dim) {
  oid_ = ++next_oid_;
  Init(n_dim);
}

Object::~Object() {
  delete position_;
  delete orientation_;
  delete velocity_;
  delete force_;
}

void Object::Init(int n_dim) {
  n_dim_ = n_dim;
  position_ = new double(n_dim_);
  orientation_ = new double(n_dim_);
  velocity_ = new double(n_dim_);
  force_ = new double(n_dim_);
  radius_ = 1;
  length_ = 0;
  ZeroForce();
}

const unsigned int Object::GetObjectID() const {
  return oid_;
}

void Object::SetRadius(double new_radius) {
  radius_ = new_radius;
}

void Object::SetLength(double new_length) {
  length_ = new_length;
}

double Object::GetRadius() const {
  return radius_;
}

double Object::GetLength() const {
  return length_;
}

void Object::ZeroForce() {
  for (int i=0; i<n_dim_; ++i)
    force_[i] = 0;
}

void Object::SetPosition(double *new_position) {
  for (int i=0; i<n_dim_; ++i)
    position_[i] = new_position[i];
}

void Object::SetOrientation(double *new_orientation) {
  for (int i=0; i<n_dim_; ++i)
    orientation_[i] = new_orientation[i];
}

void Object::SetVelocity(double *new_velocity) {
  for (int i=0; i<n_dim_; ++i)
    velocity_[i] = new_velocity[i];
}

double *Object::GetPosition() {
  return position_;
}

double *Object::GetOrientation() {
  return orientation_;
}

double *Object::GetVelocity() {
  return velocity_;
}

double *Object::GetForce() {
  return force_;
}

void Object::AddInteraction(interaction new_interaction) {
  interactions.push_back(new_interaction);
}

void Object::ApplyInteractions() {}

