#include "object.h"
#include "minimum_distance.h"

unsigned int Object::next_oid_ = 0;
unsigned int Object::next_rid_ = 0;

Object::Object(system_parameters *params, space_struct *space, long seed, SID sid) 
  : rng_(seed) {
  oid_ = ++next_oid_;
  rid_ = ++next_rid_;
  cid_ = oid_;
  sid_ = sid;
  space_ = space;
  n_dim_ = space_->n_dim;
  delta_ = params->delta;
  std::fill(position_,position_+3,0.0);
  std::fill(velocity_,velocity_+3,0.0);
  std::fill(anglevel_,anglevel_+3,0.0);
  std::fill(scaled_position_,scaled_position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  std::fill(prev_position_,prev_position_+3,0.0);
  std::fill(force_,force_+3,0.0);
  std::fill(torque_, torque_+3, 0.0);
  std::fill(dr_tot_, dr_tot_+3, 0.0);
  // Set some defaults
  diameter_ = 1;
  length_ = 0;
  k_energy_ = 0;
  p_energy_ = 0;
  is_rigid_=false;
  is_kmc_=false;
  draw_type_ = 1;
  color_ = 0;
}

int Object::DrawTypeInt(std::string dt_str) {
  if (dt_str.compare("flat") == 0)
    return 0;
  else if (dt_str.compare("orientation") == 0)
    return 1;
  else
    return 2;
}

void Object::InsertRandom() {
  double mag;
  double buffer = 0.5*diameter_;
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  // XXX Check to make sure object fits inside unit cell if non-periodic
  if (R - buffer < 0)
    error_exit("Object #%d is too large to place in system.\n system radius: %2.2f, buffer: %2.2f",GetOID(), R, buffer);
  if (space_->type.compare("sphere")==0) {
    generate_random_unit_vector(n_dim_, position_, rng_.r);
    mag = gsl_rng_uniform_pos(rng_.r) * (R - buffer);
    for (int i=0; i<n_dim_; ++i) {
      position_[i] *= mag;
    }
  }
  else if (space_->type.compare("box")==0) {
    for (int i=0; i<n_dim_; ++i)
      position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
  }
  else if (space_->type.compare("budding")==0) {
    double r = space_->bud_radius;
    double roll = gsl_rng_uniform_pos(rng_.r);
    double v_ratio;
    if (n_dim_ == 2)
      v_ratio = SQR(r) / (SQR(r)+SQR(R));
    if (n_dim_ == 3)
      v_ratio = CUBE(r) / (CUBE(r)+CUBE(R));
    mag = gsl_rng_uniform_pos(rng_.r);
    generate_random_unit_vector(n_dim_, position_, rng_.r);
    if (roll < v_ratio) {
      // Place coordinate in daughter cell
      mag *= (r - buffer);
      for (int i=0; i<n_dim_; ++i) 
        position_[i] *= mag;
      position_[n_dim_-1] += space_->bud_height;
    }
    else {
      mag *= (R - buffer);
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
    }
  }
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  UpdatePeriodic();
}

void Object::InsertRandomOriented(double *u){
  InsertRandom();
  normalize_vector(u,n_dim_);
  SetOrientation(u);
}

void Object::InsertAt(double *pos, double *u) {
  // Check to make sure position is inside unit cell if non-periodic
  double R = space_->radius;
  bool out_of_bounds = false;
  if (space_->type.compare("sphere")==0) {
    double rsq = 0;
    for (int i=0; i<n_dim_; ++i) {
      rsq += pos[i]*pos[i];
    }
    // Check that r^2 <= R^2
    if (rsq>R*R) {
      out_of_bounds = true;
    }
  }
  if (space_->type.compare("box")==0) {
    // Make sure each dimension of position, x, y etc is within the box radius
    for (int i=0; i<n_dim_; ++i) {
      if (ABS(pos[i]) > R) {
        out_of_bounds = true;
      }
    }
  }
  if (space_->type.compare("budding")==0) {
    // TODO Add checks to make sure object is placed within budding boundary...
    // right now, just trust user knows what they are doing.
  }
  if (out_of_bounds) {
    error_exit("Object %d placed outside of system unit cell! System radius: %2.2f, pos: {%2.2f %2.2f %2.2f}",GetOID(),R,pos[0],pos[1],pos[2]);
  }
  SetPosition(pos);
  normalize_vector(u,n_dim_);
  SetOrientation(u);
  UpdatePeriodic();
}

void Object::Draw(std::vector<graph_struct*> * graph_array) {
  for (int i=0;i<space_->n_periodic; ++i)
    g_.r[i] = scaled_position_[i];
  for (int i=space_->n_periodic; i<n_dim_; ++i)
    g_.r[i] = position_[i];
  std::copy(orientation_, orientation_+3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
}

// JMM - Updated to exclude updating real (global) position, since it appears not useful to do so
void Object::UpdatePeriodic() {
  double s[3], r[3];
  std::copy(scaled_position_, scaled_position_+3, s);
  periodic_boundary_conditions(space_->n_dim, space_->n_periodic, space_->unit_cell, space_->unit_cell_inv, position_, s);
  SetScaledPosition(s);
}

// Apply the force/torque/energy to this particle, override if necessary
void Object::AddForceTorqueEnergyKMC(double const * const F, double const * const T, double const p, double const k) {
    AddForce(F);
    AddTorque(T);
    AddPotential(p);
    AddKMCEnergy(k);
}

// Apply the force/torque/energy to this particle, override if necessary
void Object::AddForceTorqueEnergy(double const * const F, double const * const T, double const p) {
    AddForce(F);
    AddTorque(T);
    AddPotential(p);
}

void Object::WriteCheckpoint(std::fstream &ocheck) {
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char*>(rng_state), rng_size);
  WritePosit(ocheck);
}

void Object::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char*>(rng_state), rng_size);
  ReadPosit(icheck);
}

void Object::WritePosit(std::fstream &oposit){
  for(auto& pos : position_)
    oposit.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    oposit.write(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    oposit.write(reinterpret_cast<char*>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char*>(&length_), sizeof(length_));
}

void Object::ReadPosit(std::fstream &iposit){
  if (iposit.eof()) return;
  for(auto& pos : position_)
    iposit.read(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    iposit.read(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    iposit.read(reinterpret_cast<char*>(&u), sizeof(u));
  iposit.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char*>(&length_), sizeof(length_));
}

// Find the minimum distance beween two particles
void MinimumDistance(Simple* o1, Simple* o2, Interaction *ix, space_struct *space) {
  double const * const r1 = o1->GetRigidPosition();
  double const * const s1 = o1->GetRigidScaledPosition();
  double const * const u1 = o1->GetRigidOrientation();
  double const l1 = o1->GetRigidLength();
  double const d1 = o1->GetRigidDiameter();
  double const * const r2 = o2->GetRigidPosition();
  double const * const s2 = o2->GetRigidScaledPosition();
  double const * const u2 = o2->GetRigidOrientation();
  double const l2 = o2->GetRigidLength();
  double const d2 = o2->GetRigidDiameter();
  /* TODO: Think about how best to do this for general shapes, like 2d
     polygons that can represent the local surface of more complex 3d
     shapes. Perhaps assume all local surface to be triangular polygons.*/
  ix->dr_mag2 = 0;
  std::fill(ix->dr, ix->dr+3, 0.0);
  std::fill(ix->contact1, ix->contact1+3, 0.0);
  std::fill(ix->contact2, ix->contact2+3, 0.0);
  ix->buffer_mag = 0.5*(d1+d2);
  ix->buffer_mag2 = ix->buffer_mag*ix->buffer_mag;
  if (l1 == 0 && l2 == 0) 
    min_distance_point_point(space->n_dim, space->n_periodic, space->unit_cell,
                 r1, s1, r2, s2, ix->dr, &ix->dr_mag2);
  else if (l1 == 0 && l2 > 0)
    min_distance_sphere_sphero(space->n_dim, space->n_periodic, space->unit_cell,
                   r1, s1, r2, s2, u2, l2,
                   ix->dr, &ix->dr_mag2, ix->contact2);
  else if (l1 > 0 && l2 == 0) {
    min_distance_sphere_sphero(space->n_dim, space->n_periodic, space->unit_cell,
                   r2, s2, r1, s1, u1, l1,
                   ix->dr, &ix->dr_mag2, ix->contact1);
    for (int i=0;i<3;++i) {
      ix->dr[i] = -ix->dr[i];
    }
  }
  else if (l1 > 0 && l2 > 0) 
    min_distance_sphero(space->n_dim, space->n_periodic, space->unit_cell,
              r1, s1, u1, l1, r2, s2, u2, l2,
              ix->dr, &ix->dr_mag2, ix->contact1, ix->contact2);
}

