#include "anchor.hpp"

Anchor::Anchor() : Object() { SetSID(species_id::crosslink); }

void Anchor::Init() {
  diameter_ = params_->crosslink.diameter;
  color_ = params_->crosslink.color;
  draw_ = draw_type::_from_string(params_->crosslink.draw_type.c_str());
  Clear();
  walker_ = (params_->crosslink.walker ? true : false);
  step_direction_ = (params_->crosslink.step_direction == 0
                         ? 0
                         : SIGNOF(params_->crosslink.step_direction));
  velocity_ = params_->crosslink.velocity;
  max_velocity_ = velocity_;
  k_off_ = params_->crosslink.k_off;
  end_pausing_ = (params_->crosslink.end_pausing ? true : false);
  diffuse_ = (params_->crosslink.diffusion_flag ? true : false);
  f_spring_max_ = params_->crosslink.f_spring_max;
  f_stall_ = params_->crosslink.f_stall;
  force_dep_vel_flag_ = params_->crosslink.force_dep_vel_flag;
  SetDiffusion();
  SetSID(species_id::crosslink);
}

double const Anchor::GetMeshLambda() { return mesh_lambda_; }

double const Anchor::GetBondLambda() { return bond_lambda_; }

void Anchor::SetBondLambda(double l) { bond_lambda_ = l; }

void Anchor::SetDiffusion() { diffusion_ = sqrt(24.0 * diameter_ / delta_); }

void Anchor::SetWalker(int dir, double walk_v) {
  if (ABS(dir) != 1) {
    error_exit("Walker direction must be set to +/- 1");
  }
  walker_ = true;
  velocity_ = walk_v;
  max_velocity_ = velocity_;
  step_direction_ = dir;
}

void Anchor::UpdateAnchorPositionToMesh() {
  if (!bound_)
    return;

  // Check for unbinding of anchor from mesh
  if (gsl_rng_uniform_pos(rng_.r) <= k_off_ * delta_) {
    Clear();
    return;
  }

  // Check that the number of bonds has not changed due to dynamic instability
  if (mesh_n_bonds_ != mesh_->GetNBonds()) {
    // The number of bonds have changed, so we need to reattach to a valid bond
    bond_ = mesh_->GetBondAtLambda(mesh_lambda_);
    mesh_n_bonds_ = mesh_->GetNBonds();
  }
  /* Use the mesh to determine the bond lengths. The true bond lengths fluctuate
     about this, but should be considered approximations to the ideal mesh. */
  bond_length_ = mesh_->GetBondLength();
  mesh_length_ = mesh_->GetLength();
  /* Use current position along mesh (mesh_lambda) to determine whether the
     anchor fell off the mesh due to dynamic instability */
  if (!CheckMesh())
    return;
  // Now figure out which bond we are on in the mesh according to mesh_lambda
  bond_ = mesh_->GetBondAtLambda(mesh_lambda_);
  // Figure out how far we are from the bond tail: bond_lambda
  bond_lambda_ = mesh_lambda_ - bond_->GetBondNumber() * bond_length_;
// assert(bond_lambda_ >= 0 && bond_lambda_ <= bond_length_);
  if (bond_lambda_ < 0 || bond_lambda_ > bond_length_ + 1e-4) {
    printf("bond_num: %d\n", bond_->GetBondNumber());
    printf("mesh lambda: %2.8f\n", mesh_lambda_);
    printf("mesh length: %2.8f\n", mesh_length_);
    printf("bond lambda: %2.8f\n", bond_lambda_);
    printf("bond length: %2.8f\n", bond_length_);
    error_exit("Graargh!!!\n");
  }

  // Update anchor position with respect to bond
  UpdateAnchorPositionToBond();
}

void Anchor::UpdatePosition() {
  // Currently only bound anchors diffuse/walk (no explicit unbound anchors)
  if (!bound_ || (!diffuse_ && !walker_)) {
    return;
  }
  // Diffuse or walk along the mesh, updating mesh_lambda
  if (diffuse_) {
    Diffuse();
    if (!CheckMesh())
      return;
  }
  if (walker_) {
    Walk();
    CheckMesh();
  }
  // This occurs in UpdateAnchorPositionToMesh now
  // Now figure out which bond we are on in the mesh according to mesh_lambda
  //bond_ = mesh_->GetBondAtLambda(mesh_lambda_);
  // Figure out how far we are from the bond tail: bond_lambda
  //bond_lambda_ = mesh_lambda_ - bond_->GetBondNumber() * bond_length_;

  // assert(bond_lambda_ >= 0 && bond_lambda_ <= bond_length_);
  //if (bond_lambda_ < 0 || bond_lambda_ > bond_length_ + 1e-4) {
    //printf("bond_num: %d\n", bond_->GetBondNumber());
    //printf("mesh lambda: %2.8f\n", mesh_lambda_);
    //printf("mesh length: %2.8f\n", mesh_length_);
    //printf("bond lambda: %2.8f\n", bond_lambda_);
    //printf("bond length: %2.8f\n", bond_length_);
    //error_exit("Graargh!\n");
  //}

  // Update anchor position based on current bond attachment
  //UpdateAnchorPositionToBond();
}

void Anchor::ApplyAnchorForces() {
  if (!bound_) {
    return;
  }
  bond_->AddForce(force_);
  double dlambda[3] = {0};
  for (int i = 0; i < n_dim_; ++i) {
    dlambda[i] = (bond_lambda_ - 0.5 * bond_length_) * orientation_[i];
  }
  cross_product(dlambda, force_, torque_, 3);
  bond_->AddTorque(torque_);
}

void Anchor::Activate() {
  active_ = true;
  step_direction_ = -step_direction_;
}

void Anchor::Deactivate() {
  active_ = false;
  step_direction_ = -step_direction_;
}

void Anchor::Walk() {
  if (force_dep_vel_flag_) {
    double fmag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      fmag += force_[i] * force_[i];
    }
    fmag = sqrt(fmag);
    // Linear force-velocity relationship
    double fdep = 1 - fmag / f_stall_;
    if (fdep > 1) {
      fdep = 1;
    } else if (fdep < 0) {
      fdep = 0;
    }
    velocity_ = max_velocity_ * fdep;
  }
  double dr = step_direction_ * velocity_ * delta_;
  mesh_lambda_ += dr;
}

// Check that the anchor is still located on the filament mesh
// Returns true if anchor is still on the mesh, false otherwise
bool Anchor::CheckMesh() {
  // Check if we moved off the mesh tail
  if (mesh_lambda_ < 0) {
    // Stick to mesh ends if we have end_pausing
    if (end_pausing_) {
      mesh_lambda_ = 0;
    } else {
      // Otherwise, unbind
      Clear();
      return false;
    }
  } else if (mesh_lambda_ > mesh_length_) {
    // Same thing for walking off the head
    if (end_pausing_) {
      mesh_lambda_ = mesh_length_;
    } else {
      Clear();
      return false;
    }
  }
  return true;
}

void Anchor::Clear() {
  bound_ = false;
  bond_ = nullptr;
  mesh_ = nullptr;
  mesh_n_bonds_ = -1;
  bond_oid_ = -1;
  bond_length_ = -1;
  bond_lambda_ = -1;
  mesh_lambda_ = -1;
  active_ = false;
  ZeroForce();
}

void Anchor::Diffuse() {
  double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
  double dr = kick * diffusion_ * delta_ / diameter_;
  mesh_lambda_ += dr;
}

void Anchor::UpdateAnchorPositionToBond() {
  double const *const bond_position = bond_->GetPosition();
  double const *const bond_orientation = bond_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] = bond_orientation[i];
    position_[i] = bond_position[i] -
                   (0.5 * bond_length_ - bond_lambda_) * bond_orientation[i];
  }
  UpdatePeriodic();
}

void Anchor::Draw(std::vector<graph_struct *> *graph_array) {
  if (!bound_)
    return;
  std::copy(scaled_position_, scaled_position_ + 3, g_.r);
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_ + 3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array->push_back(&g_);
}

void Anchor::AttachObjRandom(Object *o) {
  double length = o->GetLength();
  double lambda = length * gsl_rng_uniform_pos(rng_.r);
  AttachObjLambda(o, lambda);
}

void Anchor::AttachObjLambda(Object *o, double lambda) {
  if (o->GetType() != +obj_type::bond) {
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
  bond_ = dynamic_cast<Bond *>(o);
  if (bond_ == nullptr) {
    error_exit("Object ptr passed to anchor was not referencing a bond!");
  }
  mesh_ = dynamic_cast<Mesh *>(bond_->GetMeshPtr());
  if (mesh_ == nullptr) {
    error_exit("Object ptr passed to anchor was not referencing a mesh!");
  }
  bond_oid_ = bond_->GetOID();
  mesh_n_bonds_ = mesh_->GetNBonds();
  bond_length_ = mesh_->GetBondLength();
  mesh_length_ = mesh_n_bonds_ * bond_length_;
  bond_lambda_ = lambda;

  if (bond_lambda_ < 0) {
    printf("bond_lambda: %2.2f\n", bond_lambda_);
    error_exit("Lambda passed to anchor should never be negative!");
  }
  if (bond_lambda_ > bond_length_) {
    if (bond_lambda_ - bond_length_ < 1) {
      bond_lambda_ = bond_length_;
    } else {
      error_exit(
          "Lambda passed to anchor is much larger than mesh bond length!");
    }
  }

  /* Distance anchor is relative to entire mesh length */
  mesh_lambda_ = bond_->GetBondNumber() * bond_length_ + bond_lambda_;
  SetMeshID(bond_->GetMeshID());
  UpdateAnchorPositionToBond();
  bound_ = true;
}

bool Anchor::IsBound() { return bound_; }

int const Anchor::GetBoundOID() {
  if (!bound_) {
    return -1;
  }
  return bond_->GetOID();
}

/* Temporary function for setting bound state for singly-bound crosslinks,
   in order to get them to draw while not technically bound to a bond
   ( e.g. bond_ -> null ) */
void Anchor::SetBound() { bound_ = true; }
