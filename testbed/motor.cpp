#include "motor.h"

Motor::Motor() : Site() {
  length_ = 0;
  diameter_ = 1;
  bound_ = false;
  SetDiffusion();
}

void Motor::SetDiffusion() {
  diffusion_ = sqrt(24.0*diameter_/delta_);
}

void Motor::UpdatePosition() {
  ZeroForce();
  // update motor position based on previous bond position
  //if (bound_) {
    //UpdateMotorPosition();
  //}
  // record previous position based on updated position
  //UpdatePrevPosition();
  // check probability to bind/unbind
  //if (UpdatePriors()) {
    //if (rate_kinetics_ = SLOW) {
      ////slow kinetics: if binding or unbinding, don't diffuse
      //return;
    //}
    //else if (rate_kinetics_ = FAST) {
      ////fast kinetics: diffuse even if binding or unbinding
      //Diffuse()
    //}
  //}
  //else {
    Diffuse();
  //}
  //UpdatePeriodic();
}

// update binding probabilities based on relative position to nearby bonds, then do rolls to determine if we bind or unbind
bool Motor::UpdatePriors() {
  // returns true if bound status is changed, else false
  return false;
}

void Motor::Diffuse() {
  if (bound_) {
    // If bound, diffuse along bond
    DiffuseBound();
  }
  else {
    // Otherwise diffuse normally
    for (int i=0; i<n_dim_; ++i) {
      double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
      force_[i] += kick*diffusion_;
      position_[i] += force_[i] * delta_/diameter_;
    }
  }
}

void Motor::DiffuseBound() {
  double dr[3] = {0,0,0};
  double dr_mag = 0;
  double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
  double const * const pos0 = bonds_[0].first->GetPosition();
  for (int i=0; i<n_dim_; ++i) {
    force_[i] = kick*diffusion_*orientation_[i];
  }
  for (int i = 0; i < n_dim_; ++i) {
    dr[i] = force_[i] * delta_ / diameter_;
    dr_mag += dr[i]*dr[i];
  }
  dr_mag = sqrt(dr_mag);
  // Check if we are being kicked off the edge of the bond
  bool same_bond = true;
  if (bond_lambda_ - dr_mag < 0 && kick<0) {
    //Move to previous bond if it's there
    if (same_bond = !SwitchBonds(false, dr_mag-bond_lambda_)) {
      // Otherwise move to tail of bond
      bond_lambda_ = 0.0;
    }
  }
  else if (bond_lambda_ + dr_mag > bond_length_ && kick>0) {
    //Move to next bond if it's there
    if (same_bond = !SwitchBonds(true, dr_mag - (bond_length_-bond_lambda_))) {
      // Otherwise move to head of bond
      bond_lambda_ = bond_length_;
    }
  }
  else {
    bond_lambda_ += SIGNOF(kick)*dr_mag;
  }
  // Otherwise diffuse normally along the bond
  if (same_bond) {
    for (int i=0;i<n_dim_; ++i) {
      position_[i] = pos0[i] +  (bond_lambda_ - 0.5*bond_length_)*orientation_[i];
    }
  }
}

void Motor::AttachToBond(directed_bond db, double lambda) {
  bonds_.clear();
  bonds_.push_back(std::make_pair(db.first,NONE));
  bond_length_ = db.first->GetLength();
  bond_lambda_ = (db.second == INCOMING ? bond_length_ - lambda : lambda);
  double const * const bond_position = db.first->GetPosition();
  double const * const bond_orientation = db.first->GetOrientation();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = bond_orientation[i];
    position_[i] = bond_position[i] - (0.5*bond_length_ - bond_lambda_)*orientation_[i];
  }
  bound_ = true;
}

// Returns true if switch allowed, false otherwise
bool Motor::SwitchBonds(bool next_bond, double dr_mag) {
  directed_bond db;
  if (next_bond) {
    // We went of the head, attaching to tail of next bond
    db = bonds_[0].first->GetNeighborDirectedBond(1);
  }
  else {
    // We went of the tail, attaching to head of previous bond
    db = bonds_[0].first->GetNeighborDirectedBond(0);
  }
  if (db.first==nullptr) {
    return false;
  }
  AttachToBond(db, dr_mag);
  return true;
}
