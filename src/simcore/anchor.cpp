#include "anchor.h"

Anchor::Anchor() : Object() {}

void Anchor::Init() {
  diameter_ = params_->crosslink.diameter;
  color_ = params_->crosslink.color;
  draw_ = draw_type::_from_string(params_->crosslink.draw_type.c_str());
  bound_ = false;
  walker_ = (params_->crosslink.walker ? true : false);
  step_direction_ = (params_->crosslink.step_direction == 0 ? 0 : SIGNOF(params_->crosslink.step_direction));
  velocity_ = params_->crosslink.velocity;
  max_velocity_ = velocity_;
  bond_lambda_ = 0;
  mesh_lambda_ = 0;
  active_ = false;
  diffuse_ = (params_->crosslink.diffusion_flag ? true : false);
  f_spring_max_ = params_->crosslink.f_spring_max;
  SetDiffusion();
  SetSID(species_id::crosslink);
}

double const Anchor::GetMeshLambda() {
  return mesh_lambda_;
}

void Anchor::SetDiffusion() {
  diffusion_ = sqrt(24.0*diameter_/delta_);
}

void Anchor::SetWalker(int dir, double walk_v) {
  if (ABS(dir) != 1) {
    error_exit("Walker direction must be set to +/- 1");
  }
  walker_ = true;
  velocity_ = walk_v;
  max_velocity_ = velocity_;
  step_direction_ = dir;
}

void Anchor::UpdatePosition() {
  if (!bound_) return;
  ZeroForce();
  // Update orientation based on bond if bound
  double const * const bond_orientation = bond_->GetOrientation();
  std::copy(bond_orientation, bond_orientation+3, orientation_);
  // Check for dynamic instability FIXME
  bond_length_ = bond_->GetLength();
  bond_lambda_ = mesh_lambda_ - bond_->GetBondNumber()*bond_length_ ;
  double const * const bond_position = bond_->GetPosition();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = bond_orientation[i];
    position_[i] = bond_position[i] - (0.5*bond_length_ - bond_lambda_)*orientation_[i];
  }
  if (walker_) {
    Walk();
  }
  //if (!anchored_) {
    //// This can change anchored_
    //CheckNearBoundary();
  //}
  //if (anchored_) {
    //ApplyAnchorForces();
  //}
  if (diffuse_) {
    Diffuse();
  }
  /* Clear forces, since they will be used for storing crosslink tether 
   * forces */
  ZeroForce();
  UpdatePeriodic();
}

void Anchor::ApplyAnchorForces() {
  /* TODO add torques */
  bond_->AddForce(force_);
}

//void Anchor::ApplyAnchorForces() {
  //double dr[3] = {0, 0, 0};
  //double temp[3] = {0, 0, 0};
  //double f_mag = 0.0;
  //for (int i=0;i<n_dim_; ++i) {
    //dr[i] = position_[i] - anchor_.position_[i];
    //anchor_.force_[i] = - anchor_.k_spring_ * dr[i];
    //temp[i] = (bond_lambda_ - 0.5*bond_length_)*orientation_[i];
    //f_mag += anchor_.force_[i]*anchor_.force_[i];
  //}
  //if (f_mag > SQR(f_spring_max_)) {
    //DetachBoundary();
    //return;
  //}
  //f_mag = sqrt(f_mag);
  //cross_product(temp, anchor_.force_, anchor_.torque_, n_dim_);
  //bonds_[0].first->AddForce(anchor_.force_);
  //bonds_[0].first->AddTorque(anchor_.torque_);
  //velocity_ = max_velocity_ * exp(-pow(f_mag/f_spring_max_, 4));
//}

//void Anchor::CheckNearBoundary() {
  //// Only attaches to spherical and budding yeast boundary conditions
  //if (params_->boundary == 3) {
    //CheckNearBuddingBoundary();
  //}
  //else if (params_->boundary == 2) {
    //double dr_mag = 0.0;
    //for (int i=0;i<n_dim_;++i) {
      //dr_mag += position_[i]*position_[i];
    //}
    //if (dr_mag > SQR(space_->radius - 0.5*diameter_ - 1)) {
      //AnchorBoundary(position_);
    //}
  //}
//}

// Copy of MinimumDistance::PointBuddingBC (see for comments/notes)
//void Anchor::CheckNearBuddingBoundary() {
  //double const * const r = position_;
  //double dr[3] = {0, 0, 0};
  //double dr_mag2 = 0;
  //bool in_mother = (r[n_dim_-1] < space_->bud_neck_height);
  //bool in_cone_region = !( r[n_dim_-1] < 0 || r[n_dim_-1] > space_->bud_height );
  //double r_mag = 0.0;
  //for (int i=0;i<n_dim_-1;++i) {
    //r_mag += SQR(r[i]);
  //}
  //double z0 = (in_mother ? 0 : space_->bud_height);
  //[> Check if in_cone_region means what it says it means <]
  //if (in_cone_region) {
    //double cone_rho2 = SQR(space_->bud_neck_radius)*SQR(r[n_dim_-1]-z0)/SQR(space_->bud_neck_height);
    //in_cone_region = (r_mag < cone_rho2);
  //}
  //[> Now in_cone_region definitely means what it says it means <]
  //if (in_cone_region) {
    //double scale_factor = space_->bud_neck_radius/sqrt(r_mag) - 1;
    //if (scale_factor < 0) error_exit("Something went wrong in CheckNearBuddingBoundary !\n");
    //dr_mag2 = 0;
    //for (int i=0;i<n_dim_-1; ++i) {
      //dr[i] = scale_factor*r[i];
      //dr_mag2 += SQR(dr[i]);
    //}
    //dr[n_dim_-1] = space_->bud_neck_height - r[n_dim_-1];
    //dr_mag2 += SQR(dr[n_dim_-1]);
  //}
  //else {
    //r_mag = sqrt(r_mag+SQR(r[n_dim_-1] - z0));
    //double r_cell = (in_mother ? space_->radius : space_->bud_radius);
    //dr_mag2 = 0;
    //for (int i=0;i<n_dim_-1;++i) {
      //dr[i] = ((r_cell-diameter_)/r_mag-1)*r[i];
      //dr_mag2 += SQR(dr[i]);
    //}
    //dr[n_dim_-1] = ((r_cell-diameter_)/r_mag-1)*(r[n_dim_-1] - z0);
    //dr_mag2 += SQR(dr[n_dim_-1]);
  //}
  //if (dr_mag2 < SQR(1 + 0.5*diameter_)) {
    //AnchorBoundary(position_);
  //}
//}


//void Anchor::AnchorBoundary(double * attach_point) {
  //anchored_ = true;
  //Activate();
  //for (int i=0;i<n_dim_; ++i) {
    //anchor_.position_[i] = position_[i];
  //}
//}

//void Anchor::DetachBoundary() {
  //anchored_ = false;
  //velocity_ = max_velocity_;
  //Deactivate();
//}


void Anchor::Activate() {
  active_ = true;
  step_direction_ = -step_direction_;
}

void Anchor::Deactivate() {
  active_ = false;
  step_direction_ = -step_direction_;
}

void Anchor::Walk() {
  double dr[3] = {0, 0, 0};
  double dr_mag = velocity_*delta_;
  double const * const pos0 = bond_->GetPosition();
  for (int i = 0; i < n_dim_; ++i) {
    dr[i] = step_direction_ * dr_mag * orientation_[i];
  }
  // Check if we are walking off the edge of the bond
  bool same_bond = true;
  if (bond_lambda_ - dr_mag < 0 && step_direction_<0) {
    //Move to previous bond if it's there
    if (same_bond = !SwitchBonds(false, dr_mag-bond_lambda_)) {
      // Otherwise move to tail of bond
      mesh_lambda_ -= bond_lambda_;
      bond_lambda_ = 0.0;
    }
  }
  else if (bond_lambda_ + dr_mag > bond_length_ && step_direction_>0) {
    //Move to next bond if it's there
    if (same_bond = !SwitchBonds(true, dr_mag - (bond_length_-bond_lambda_))) {
      // Otherwise move to head of bond
      mesh_lambda_ += bond_length_ - bond_lambda_;
      bond_lambda_ = bond_length_;
    }
  }
  else {
    bond_lambda_ += step_direction_*dr_mag;
    mesh_lambda_ += step_direction_*dr_mag;
  }
  // Otherwise step normally along the bond
  if (same_bond) {
    for (int i=0;i<n_dim_; ++i) {
      position_[i] = pos0[i] +  (bond_lambda_ - 0.5*bond_length_)*orientation_[i];
    }
  }
  else {
    mesh_lambda_ += step_direction_*dr_mag;
  }
}

//void Anchor::Diffuse() {
  //if (bound_) {
    //// If bound, diffuse along bond
    //DiffuseBound();
  //}
  //else {
    //// Otherwise diffuse normally
    //for (int i=0; i<n_dim_; ++i) {
      //double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
      //force_[i] += kick*diffusion_;
      //position_[i] += force_[i] * delta_/diameter_;
    //}
  //}
//}

// Formerly DiffuseBound

void Anchor::Clear() {
  bound_ = false;
  bond_ = nullptr;
  bond_length_ = -1;
  bond_lambda_ = -1;
  mesh_lambda_ = -1;
  ZeroForce();
}

void Anchor::Diffuse() {
  double dr[3] = {0, 0, 0};
  double dr_mag = 0;
  double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
  double const * const pos0 = bond_->GetPosition();
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
      mesh_lambda_ -= bond_lambda_;
      bond_lambda_ = 0.0;
    }
  }
  else if (kick > 0 && bond_lambda_ + dr_mag > (bond_length_ = bond_->GetLength()) && kick>0) {
    //Move to next bond if it's there
    if (same_bond = !SwitchBonds(true, dr_mag - (bond_length_-bond_lambda_))) {
      // Otherwise move to head of bond
      mesh_lambda_ += bond_length_ - bond_lambda_;
      bond_lambda_ = bond_length_;
    }
  }
  else {
    double dr = SIGNOF(kick)*dr_mag;
    bond_lambda_ += dr;
    mesh_lambda_ += dr;
  }
  // Otherwise diffuse normally along the bond
  if (same_bond) {
    for (int i=0;i<n_dim_; ++i) {
      position_[i] = pos0[i] +  (bond_lambda_ - 0.5*bond_length_)*orientation_[i];
    }
  }
  else {
    mesh_lambda_ += SIGNOF(kick)*dr_mag;
  }
}

void Anchor::AttachToBond(directed_bond db, double lambda, double mesh_lambda) {
  bond_ = db.first;
  bond_length_ = db.first->GetLength();
  bond_lambda_ = (db.second == INCOMING ? bond_length_ - lambda : lambda);
  mesh_lambda_ = mesh_lambda;
  double const * const bond_position = db.first->GetPosition();
  double const * const bond_orientation = db.first->GetOrientation();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = bond_orientation[i];
    position_[i] = bond_position[i] - (0.5*bond_length_ - bond_lambda_)*orientation_[i];
  }
  bound_ = true;
  UpdatePeriodic();
}

void Anchor::Draw(std::vector<graph_struct*> * graph_array) {
  if (!bound_) return;
  std::copy(scaled_position_, scaled_position_+3, g_.r);
  for (int i=space_->n_periodic; i<n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_+3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array->push_back(&g_);
  //printf("%2.2f %2.2f %2.2f\n", g_.r[0], g_.r[1], g_.r[2]);
  //printf("%2.2f %2.2f %2.2f\n", g_.u[0], g_.u[1], g_.u[2]);
  //printf("%2.2f %2.2f\n\n", g_.length, g_.diameter);
}



// Returns true if switch allowed, false otherwise
bool Anchor::SwitchBonds(bool next_bond, double dr_mag) {
  directed_bond db;
  if (next_bond) {
    // We went off the head, attaching to tail of next bond
    db = bond_->GetNeighborDirectedBond(1);
  }
  else {
    // We went off the tail, attaching to head of previous bond
    db = bond_->GetNeighborDirectedBond(0);
  }
  if (db.first==nullptr) {
    return false;
  }
  AttachToBond(db, dr_mag, mesh_lambda_);
  return true;
}

void Anchor::AttachObjRandom(Object * o) {
  if (o->GetType() != +obj_type::bond) {
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
  bond_ = dynamic_cast<Bond*>(o);
  bond_length_ = bond_->GetLength();
  bond_lambda_ = bond_length_*gsl_rng_uniform_pos(rng_.r);
  /* First bond has bond number of zero */
  int bond_number = bond_->GetBondNumber();
  /* Distance anchor is relative to entire mesh length */
  mesh_lambda_ = bond_number*bond_length_ + bond_lambda_;
  directed_bond db = std::make_pair(bond_, OUTGOING);
  /* Passing these is redundant, but perhaps it is important for
     AttachToBond to take these arguments as a public method? */
  AttachToBond(db, bond_lambda_, mesh_lambda_);
  SetMeshID(bond_->GetMeshID());
}

bool Anchor::IsBound() {
  return bound_;
}
