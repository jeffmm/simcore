#include "object.h"
#include "output_manager.h"

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
  interactions_.clear();
  // Set some defaults
  diameter_ = 1;
  length_ = 0;
  k_energy_ = 0;
  p_energy_ = 0;
  is_rigid_=false;
  is_kmc_=false;
  draw_type_ = 1;
  color_[0] = 1.0;
  color_[1] = 0.0;
  color_[2] = 0.0;
  color_[3] = 1.0;
}

//void Object::Dump() {
  //printf("  oid %d:\n",oid_);
  //printf("    sid %d\n    cid %d\n    rid %d\n",sid_,cid_,rid_);
  //printf("    r : {%f, %f, %f}\n", position_[0], position_[1], position_[2]);
  //printf("    s : {%f, %f, %f}\n", scaled_position_[0], scaled_position_[1], scaled_position_[2]);
  //printf("    u : {%f, %f, %f}\n", orientation_[0], orientation_[1], orientation_[2]);
  //printf("    vel : {%f, %f, %f}\n", velocity_[0], velocity_[1], velocity_[2]);
  //printf("    prev r : {%f, %f, %f}\n", prev_position_[0], prev_position_[1], prev_position_[2]);
  //printf("    f : {%f, %f, %f}\n", force_[0], force_[1], force_[2]);
  //printf("    t : {%f, %f, %f}\n", torque_[0], torque_[1], torque_[2]);
  //printf("    d: %f\n    l: %f\n    ke: %f\n    pe: %f\n",diameter_,length_,k_energy_,p_energy_);
//}

void Object::InsertRandom(double buffer) {
  double mag;
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  if (R - buffer < 0)
    error_exit("ERROR: Object #%d is too large to place in system.\n system radius: %2.2f, buffer: %2.2f\n",GetOID(), R, buffer);
  if (space_->type.compare("sphere")==0) {
    generate_random_unit_vector(n_dim_, position_, rng_.r);
    mag = gsl_rng_uniform_pos(rng_.r) * (R - buffer);
    for (int i=0; i<n_dim_; ++i) {
      position_[i] *= mag;
    }
  }
  else if (space_->type.compare("cube")==0) {
    for (int i=0; i<n_dim_; ++i)
      position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer);
  }
  else if (space_->type.compare("snowman")==0) {
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
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
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

void Object::InsertOriented(double* buffer, const double* const u){
  SetOrientation(u);
  double mag;
  if (space_->n_periodic == n_dim_)
    std::fill(buffer, buffer+3 , 0);
  double R = space_->radius;
  double buffer_mag = sqrt(dot_product(n_dim_, buffer, buffer));
  if (R - buffer_mag < 0) 
    error_exit("ERROR: Object #%d is too large to place in system.\n",GetOID());
  if (space_->type.compare("sphere")==0) {
    generate_random_unit_vector(n_dim_, position_, rng_.r);
    mag = gsl_rng_uniform_pos(rng_.r) * (R - buffer_mag);
    for (int i=0; i<n_dim_; ++i) {
      position_[i] *= mag;
    }
  }
  else if (space_->type.compare("cube")==0) {
    for (int i=0; i<n_dim_; ++i)
      position_[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (R - buffer[i]);
  }
  else if (space_->type.compare("snowman")==0) {
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
      mag *= (r - buffer_mag);
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
      position_[n_dim_-1] += space_->bud_height;
    }
    else {
      mag *= (R - buffer_mag);
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
      }
    }
  }
  UpdatePeriodic();
}

void Object::Draw(std::vector<graph_struct*> * graph_array) {
  std::copy(position_, position_+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  std::copy(color_, color_+4, g_.color);
  //g_.length = (length_-diameter_ > 0 ? length_-diameter_ : 0);
  g_.length = length_;
  g_.diameter = diameter_;
  g_.draw_type = draw_type_;
  graph_array->push_back(&g_);
}

void Object::UpdatePeriodic() {
  if (space_->n_periodic == 0)
    return;
  double r[3], s[3], rp[3], drp[3];
  std::copy(orientation_, orientation_+3, g_.u);
  std::copy(position_, position_+3, r);
  std::copy(scaled_position_, scaled_position_+3, s);
  std::copy(prev_position_, prev_position_+3, rp);
  for (int i=0; i<n_dim_; ++i)
    drp[i] = r[i]-rp[i];
  periodic_boundary_conditions(space_->n_periodic, space_->unit_cell, space_->unit_cell_inv, r, s);
  SetPosition(r);
  SetScaledPosition(s);
  for (int i=0; i<n_dim_; ++i)
    rp[i] = r[i]-drp[i];
  SetPrevPosition(rp);
}

void Simple::ApplyInteractions() {
  for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
    auto ix=(*it);
    double f[4];
    //ix.potential->CalcPotential(ix.dr, ix.dr_mag, ix.buffer, f);
    AddForce(f);
    AddPotential(f[n_dim_]);
  }
}

// Apply the force/torque/energy to this particle, override if necessary
void Object::AddForceTorqueEnergyKMC(double const * const F, double const * const T, double const p, double const k) {
    // XXX: CJE assume that we need to zero out the force and energy first
    //ZeroForce();
    AddForce(F);
    AddTorque(T);
    AddPotential(p);
    AddKMCEnergy(k);
}

// Apply the force/torque/energy to this particle, override if necessary
void Object::AddForceTorqueEnergy(double const * const F, double const * const T, double const p) {
    // XXX: CJE assume that we need to zero out the force and energy first
    //ZeroForce();
    AddForce(F);
    AddTorque(T);
    AddPotential(p);
}

// Find the minimum distance beween two particles
void MinimumDistance(Simple* o1, Simple* o2, interactionmindist& imd, int& ndim, int& nperiodic, space_struct *space) {
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
  imd.dr_mag2 = 0;
  std::fill(imd.dr, imd.dr+3, 0.0);
  std::fill(imd.contact1, imd.contact1+3, 0.0);
  std::fill(imd.contact2, imd.contact2+3, 0.0);
  imd.buffer_mag = 0.5*(d1+d2);
  imd.buffer_mag2 = imd.buffer_mag*imd.buffer_mag;
  if (l1 == 0 && l2 == 0) 
    min_distance_point_point(ndim, nperiodic, space->unit_cell,
                 r1, s1, r2, s2, imd.dr, &imd.dr_mag2);
  else if (l1 == 0 && l2 > 0)
    min_distance_sphere_sphero(ndim, nperiodic, space->unit_cell,
                   r1, s1, r2, s2, u2, l2,
                   imd.dr, &imd.dr_mag2, imd.contact2);
  else if (l1 > 0 && l2 == 0)
    min_distance_sphere_sphero(ndim, nperiodic, space->unit_cell,
                   r2, s2, r1, s1, u1, l1,
                   imd.dr, &imd.dr_mag2, imd.contact1);
  else if (l1 > 0 && l2 > 0) 
    min_distance_sphero(ndim, nperiodic, space->unit_cell,
              r1, s1, u1, l1, r2, s2, u2, l2,
              imd.dr, &imd.dr_mag2, imd.contact1, imd.contact2);
  //imd.dr_mag = sqrt(imd.dr_mag2); // XXX FIXME should take out and only compute if necessary
}

void Object::WritePosit(std::fstream &op){
  for(auto& pos : position_)
    op.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    op.write(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    op.write(reinterpret_cast<char*>(&u), sizeof(u));
  op.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  op.write(reinterpret_cast<char*>(&length_), sizeof(length_));
}

void Object::ReadPosit(std::fstream &ip){
  if (ip.eof()) return;
  for(auto& pos : position_)
    ip.read(reinterpret_cast<char*>(&pos), sizeof(pos));
  for(auto& spos : scaled_position_)
    ip.read(reinterpret_cast<char*>(&spos), sizeof(spos));
  for(auto& u : orientation_)
    ip.read(reinterpret_cast<char*>(&u), sizeof(u));
  ip.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ip.read(reinterpret_cast<char*>(&length_), sizeof(length_));
}

//void Object::AddVirial(double const * const f, double const * const dr) {
//void Object::AddVirial(double *f, double *dr) {
  //for (int i=0; i<n_dim_; ++i)
    //for (int j=i; j<n_dim_; ++j){
      //printf("    fr%d=%f, dr%d=%f\n",i,f[i],j,dr[j]);
      //spec_virial_[3*i+j] = spec_virial_[3*j+i] += f[i] * dr[j];
    //}
//}



