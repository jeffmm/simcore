#include "object.h"

/****************************
** Object member functions **
****************************/
Object::Object(long seed) : rng_(seed) {
  oid_ = _next_oid_++; 
  if (params_ == nullptr) {
    std::cerr << "ERROR! Object parameters not set before instantiation!\n";
    exit(1);
  }
  if (n_dim_ == 0) {
    std::cerr << "ERROR! Object dimensions not set before instantiation!\n";
    exit(1);
  }
  if (delta_ == 0) {
    std::cerr << "ERROR! Object delta not set before instantiation!\n";
    exit(1);
  }
  std::fill(position_,position_+3,0.0);
  std::fill(scaled_position_,scaled_position_+3,0.0);
  std::fill(prev_position_,prev_position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  std::fill(force_,force_+3,0.0);
  std::fill(torque_,torque_+3,0.0);
  length_ = 0;
  diameter_ = 0;
}

int Object::_next_oid_ = 0;
long Object::_seed_ = 7777777;
int Object::n_dim_ = 0;
double Object::delta_ = 0;
system_parameters const * Object::params_ = nullptr;
int Object::GetOID() {
  return oid_;
}
double const * const Object::GetPosition() {
  return position_;
}
double const * const Object::GetOrientation() {
  return orientation_;
}
double const Object::GetLength() {
  return length_;
}
double const Object::GetDiameter() {
  return diameter_;
}
void Object::SetPosition(double * pos) {
  std::copy (pos, pos+3, position_);
}
void Object::SetOrientation(double * u) {
  std::copy (u, u+3, orientation_);
}
void Object::SetLength(double l) {
  length_ = l;
}
void Object::SetDiameter(double d) {
  diameter_ = d;
}
void Object::ZeroForce() {
  std::fill(force_,force_+3,0.0);
  std::fill(torque_,torque_+3,0.0);
}
void Object::UpdatePosition() {
  UpdatePrevPosition();
  // update object position;
}
void Object::UpdatePrevPosition() {
  for (int i=0; i<n_dim_; ++i) {
    prev_position_[i] = position_[i];
  }
}
void Object::Report() {
  fprintf(stderr, "    OID: %d\n",GetOID());
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n",position_[0],position_[1],position_[2]);
  fprintf(stderr, "      u: {%2.2f %2.2f %2.2f}\n",orientation_[0],orientation_[1],orientation_[2]);
  fprintf(stderr, "      l: %2.2f\n",length_);
  fprintf(stderr, "      d: %2.2f\n",diameter_);
}
// Main draw function, return struct of graphics info
void Object::Draw(std::vector<graph_struct*> * graph_array) {
  // Record current position, etc in graph_struct
  for (int i=0; i<3; ++i) {
    g_struct.r[i] = position_[i];
    g_struct.u[i] = orientation_[i];
  }
  g_struct.length = length_;
  g_struct.diameter = diameter_;
  g_struct.draw_type = 1;
  // push graph_struct into array
  graph_array->push_back(&g_struct);
}

