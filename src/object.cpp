#include "object.h"

unsigned int Object::next_oid_ = 0;
unsigned int Object::next_cid_ = 0;

Object::Object(system_parameters *params, space_struct *space, long seed, SID sid) {
  oid_ = ++next_oid_;
  sid_ = sid;
  space_ = space;
  n_dim_ = space_->n_dim; 
  delta_ = params->delta;
  std::fill(position_,position_+3,0.0);
  std::fill(scaled_position_,scaled_position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  std::fill(prev_position_,prev_position_+3,0.0);
  std::fill(force_,force_+3,0.0);
  rng_.init(seed);
  interactions_.clear();
  // Set some defaults
  diameter_ = 1;
  length_ = 0;
  k_energy_ = 0;
  is_simple_=false;
}

Object::Object(const Object& that) {
  n_dim_=that.n_dim_;
  oid_ = that.oid_;
  cid_ = that.cid_;
  sid_ = that.sid_;
  next_oid_=that.next_oid_;
  delta_ = that.delta_;
  std::copy(that.position_, that.position_+3, position_);
  std::copy(that.scaled_position_, that.scaled_position_+3, scaled_position_);
  std::copy(that.prev_position_, that.prev_position_+3, prev_position_);
  std::copy(that.orientation_, that.orientation_+3, orientation_);
  std::copy(that.force_, that.force_+3, force_);
  diameter_ = that.diameter_;
  length_ = that.length_;
  k_energy_ = that.k_energy_;
  p_energy_ = that.p_energy_;
  space_ = that.space_;
  interactions_ = that.interactions_;
  g_ = that.g_;
  rng_.init(gsl_rng_get(that.rng_.r));
  is_simple_=that.is_simple_;
}

Object &Object::operator=(Object const& that) {
  n_dim_=that.n_dim_;
  oid_=that.oid_;
  cid_=that.cid_;
  sid_=that.sid_;
  next_oid_=that.next_oid_;
  delta_ = that.delta_;
  std::copy(that.position_, that.position_+3, position_);
  std::copy(that.scaled_position_, that.scaled_position_+3, scaled_position_);
  std::copy(that.prev_position_, that.prev_position_+3, prev_position_);
  std::copy(that.orientation_, that.orientation_+3, orientation_);
  std::copy(that.force_, that.force_+3, force_);
  diameter_ = that.diameter_;
  length_ = that.length_;
  k_energy_ = that.k_energy_;
  p_energy_ = that.p_energy_;
  space_ = that.space_;
  interactions_ = that.interactions_;
  g_ = that.g_;
  rng_.init(gsl_rng_get(that.rng_.r));
  is_simple_=that.is_simple_;
  return *this;
}

void Object::InsertRandom(double buffer) {
  double mag;
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
      mag *= (space_->radius - buffer);
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
    #ifdef DEBUG
    printf("o(%d) = ", GetOID());
    printf("{%2.2f, %2.2f}, ", GetPosition()[0], GetPosition()[1]);
    printf("f{%2.2f, %2.2f, %2.2f}\n", force_[0], force_[1], p_energy_);
    #endif
  }
}
