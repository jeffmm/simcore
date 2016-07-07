#include "object.h"

#include "minimum_distance.h"

unsigned int Object::next_oid_ = 0;
unsigned int Object::next_rid_ = 0;

Object::Object(system_parameters *params, space_struct *space, long seed, SID sid) {
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
  rng_.init(seed);
  interactions_.clear();
  // Set some defaults
  diameter_ = 1;
  length_ = 0;
  k_energy_ = 0;
  p_energy_ = 0;
  is_rigid_=false;
  is_kmc_=false;
}

Object::Object(const Object& that) {
  n_dim_=that.n_dim_;
  oid_ = that.oid_;
  cid_ = that.cid_;
  sid_ = that.sid_;
  rid_ = that.rid_;
  next_oid_=that.next_oid_;
  delta_ = that.delta_;
  std::copy(that.position_, that.position_+3, position_);
  std::copy(that.velocity_, that.velocity_+3, velocity_);
  std::copy(that.anglevel_, that.anglevel_+3, anglevel_);
  std::copy(that.scaled_position_, that.scaled_position_+3, scaled_position_);
  std::copy(that.prev_position_, that.prev_position_+3, prev_position_);
  std::copy(that.orientation_, that.orientation_+3, orientation_);
  std::copy(that.force_, that.force_+3, force_);
  std::copy(that.torque_, that.torque_+3, torque_);
  std::copy(that.dr_tot_, that.dr_tot_+3, dr_tot_);
  diameter_ = that.diameter_;
  length_ = that.length_;
  k_energy_ = that.k_energy_;
  p_energy_ = that.p_energy_;
  space_ = that.space_;
  interactions_ = that.interactions_;
  g_ = that.g_;
  rng_.init(gsl_rng_get(that.rng_.r));
  is_rigid_=that.is_rigid_;
  is_kmc_=that.is_kmc_;
}

Object &Object::operator=(Object const& that) {
  n_dim_=that.n_dim_;
  oid_=that.oid_;
  cid_=that.cid_;
  sid_=that.sid_;
  rid_=that.rid_;
  next_oid_=that.next_oid_;
  delta_ = that.delta_;
  std::copy(that.position_, that.position_+3, position_);
  std::copy(that.velocity_, that.velocity_+3, velocity_);
  std::copy(that.anglevel_, that.anglevel_+3, anglevel_);
  std::copy(that.scaled_position_, that.scaled_position_+3, scaled_position_);
  std::copy(that.prev_position_, that.prev_position_+3, prev_position_);
  std::copy(that.orientation_, that.orientation_+3, orientation_);
  std::copy(that.force_, that.force_+3, force_);
  std::copy(that.torque_, that.torque_+3, torque_);
  std::copy(that.dr_tot_, that.dr_tot_+3, dr_tot_);
  diameter_ = that.diameter_;
  length_ = that.length_;
  k_energy_ = that.k_energy_;
  p_energy_ = that.p_energy_;
  space_ = that.space_;
  interactions_ = that.interactions_;
  g_ = that.g_;
  rng_.init(gsl_rng_get(that.rng_.r));
  is_rigid_=that.is_rigid_;
  is_kmc_=that.is_kmc_;
  return *this;
}

void Object::InsertRandom(double buffer) {
  double mag;
  if (space_->n_periodic == n_dim_)
    buffer = 0;
  double R = space_->radius;
  if (R - buffer < 0) 
    error_exit("ERROR: Object #%d is too large to place in system.\n",GetOID());
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

void Object::Draw(std::vector<graph_struct*> * graph_array) {
  std::copy(position_, position_+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  //g_.length = (length_-diameter_ > 0 ? length_-diameter_ : 0);
  g_.length = length_;
  g_.diameter = diameter_;
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
    ix.potential->CalcPotential(ix.dr, ix.dr_mag, ix.buffer, f);
    AddForce(f);
    AddPotential(f[n_dim_]);
  }
}

// Apply the force/torque/energy to this particle, override if necessary
void Object::AddForceTorqueEnergyKMC(double const * const F, double const * const T, double const p, double const k) {
    // XXX: CJE assume that we need to zero out the force and energy first
    ZeroForce();
    AddForce(F);
    AddTorque(T);
    AddPotential(p);
    AddKMCEnergy(k);
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
  imd.dr_mag = sqrt(imd.dr_mag2);
}
