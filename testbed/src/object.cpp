#include "object.h"
void generate_random_unit_vector(double * vec, int const n_dim, gsl_rng * r);
void get_random_coordinate(double * pos, int const n_dim, double radius, boundary_type btype, gsl_rng * r);


/****************************
** Object member functions **
****************************/
Object::Object() {
  oid_ = _next_oid_++; 
  if (_params_ == nullptr) {
    std::cerr << "ERROR! Object parameters not set before instantiation!\n";
    exit(1);
  }
  params_ = _params_;
  if (_n_dim_ == 0) {
    std::cerr << "ERROR! Object dimensions not set before instantiation!\n";
    exit(1);
  }
  n_dim_ = _n_dim_;
  std::fill(position_,position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  length_ = 0;
  diameter_ = 0;
  rng_.Init(_seed_);
  _seed_ = gsl_rng_get(rng_.r);
}

int Object::_next_oid_ = 0;
int Object::_n_dim_ = 0;
long Object::_seed_ = 7777777;
system_parameters const * Object::_params_ = nullptr;
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

/**************************
** Site member functions **
**************************/
void Site::AddBond(Bond * bond) {
  bonds_.push_back(bond);
}
Bond * Site::GetBond(int i) {
  if (i<0 || i >= bonds_.size()) {
    std::cerr << "ERROR! Requested adjacent bond out of bounds!\n";
  }
  return bonds_[i];
}
Bond * Site::GetOtherBond(int bond_oid) {
  for (std::vector<Bond*>::iterator it=bonds_.begin(); it!=bonds_.end(); ++it) {
    if ((*it)->GetOID() != bond_oid) {
      return (*it);
    }
  }
  // There are no other bonds
  return nullptr;
}
void Site::Report() {
  fprintf(stderr, "  Site:\n");
  fprintf(stderr, "    OID: %d\n",GetOID());
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n",position_[0],position_[1],position_[2]);
  fprintf(stderr, "      d: %2.2f\n",diameter_);
}

/**************************
** Bond member functions **
**************************/
Bond::Bond(Site * s1, Site * s2) {
  Init(s1,s2);
}
void Bond::Init(Site * s1, Site * s2) {
  s1->AddBond(this);
  s2->AddBond(this);
  sites_[0] = s1;
  sites_[1] = s2;
  double const * const r1 = s1->GetPosition();
  double const * const r2 = s2->GetPosition();
  diameter_ = s1->GetDiameter();
  length_ = 0;
  for (int i=0;i<n_dim_;++i) {
    orientation_[i] = r2[i] - r1[i];
    length_ += orientation_[i]*orientation_[i];
  }
  length_ = sqrt(length_);
  for (int i=0;i<n_dim_;++i) {
    position_[i] = r1[i] + 0.5*orientation_[i];
    orientation_[i]/=length_;
  }
}
Site * Bond::GetSite(int i) {
  if (i<0 || i>1) {
    std::cerr << "ERROR! Requested adjacent site out of bounds!\n";
  }
  return sites_[i];
}
Bond * Bond::GetNeighborBond(int i) {
  if (i<0 || i>1) {
    std::cerr << "ERROR! Requested neighboring bond out of bounds!\n";
  }
  return sites_[i]->GetOtherBond(GetOID());
}
void Bond::Report() {
  fprintf(stderr, "  Bond:\n");
  Object::Report();
}

/**************************
** Mesh member functions **
**************************/
void Mesh::AddSite(Site s) {
  sites_.push_back(s);
  n_sites_++;
}
void Mesh::AddBond(Bond b) {
  bonds_.push_back(b);
  n_bonds_++;
}

void Mesh::InitSiteAt(double * pos, double d) {
  Site s;
  s.SetPosition(pos);
  s.SetDiameter(d);
  AddSite(s);
}
void Mesh::InitBondAt(double * pos, double * u, double l, double d) {
  Site s1,s2;
  s1.SetDiameter(d);
  s2.SetDiameter(d);
  for (int i=0; i<n_dim_; ++i) {
    pos[i] -= 0.5*l*u[i];
  }
  s1.SetPosition(pos);
  for (int i=0; i<n_dim_; ++i) {
    pos[i] += l*u[i];
  }
  s2.SetPosition(pos);
  AddSite(s1);
  AddSite(s2);
  Bond b(&sites_[n_sites_-2],&sites_[n_sites_-1]);
  AddBond(b);
}
void Mesh::InitRandomSite(double d) {
  double pos[3];
  get_random_coordinate(pos,params_->n_dim,params_->system_radius,params_->boundary,rng_.r);
  InitSiteAt(pos,d);
}
// Default d=1
void Mesh::AddRandomBondAnywhere(double l, double d) {
  if (n_sites_ < 1) {
    InitRandomSite(d);
  }
  int i_site=gsl_rng_uniform_int(rng_.r,n_sites_);
  AddRandomBondToSite(l, i_site);
}
void Mesh::AddRandomBondToSite(double l, int i_site) {
  if (i_site > n_sites_ || i_site < 0) {
    std::cerr << "ERROR! Site index out of range in AddRandomBondToSite!\n";
    exit(1);
  }
  double const d = sites_[i_site].GetDiameter();
  double const * const pos0 = sites_[i_site].GetPosition();
  double pos[3];
  generate_random_unit_vector(pos, n_dim_, rng_.r);
  for (int i=0; i<n_dim_; ++i) {
    pos[i] = pos0[i] + l * pos[i];
  }
  InitSiteAt(pos,d);
  Bond b(&sites_[i_site],&sites_[n_sites_-1]);
  AddBond(b);
}
void Mesh::AddRandomBondToTip(double l) {
  AddRandomBondToSite(l,n_sites_-1);
}
void Mesh::AddBondToTip(double *u, double l) {
  AddBondToSite(u,l,n_sites_-1);
}
void Mesh::AddBondToSite(double *u, double l, int i_site) {
  if (i_site > n_sites_ || i_site < 0) {
    std::cerr << "ERROR! Site index out of range in AddBondToSite!\n";
    exit(1);
  }
  double const d = sites_[i_site].GetDiameter();
  double const * const pos0 = sites_[i_site].GetPosition();
  double pos[3];
  for (int i=0; i<n_dim_; ++i) {
    pos[i] = pos0[i] + l * u[i];
  }
  InitSiteAt(pos,d);
  Bond b(&sites_[i_site],&sites_[n_sites_-1]);
  AddBond(b);
}
void Mesh::ReportSites() {
  for (std::vector<Site>::iterator it=sites_.begin(); it!=sites_.end(); ++it) {
    it->Report();
  }
}
void Mesh::ReportBonds() {
  for (std::vector<Bond>::iterator it=bonds_.begin(); it!=bonds_.end(); ++it) {
    it->Report();
  }
}
void Mesh::Report() {
  fprintf(stderr,"Mesh: \n");
  fprintf(stderr,"  OID: %d\n",GetOID());
  fprintf(stderr,"  n_sites: %d\n",n_sites_);
  fprintf(stderr,"  n_bonds: %d\n",n_bonds_);
  ReportSites();
  ReportBonds();
}
void Mesh::SetBondLength(double l) {
  bond_length_ = l;
}
void Mesh::Draw(std::vector<graph_struct*> * graph_array) {
  for (std::vector<Bond>::iterator it = bonds_.begin(); it!=bonds_.end(); ++it) {
    it->Draw(graph_array);
  }
}

Motor::Motor() {
  bound_ = false;
}

void Motor::UpdatePosition() {
  // update motor position based on previous bond position
  if (bound_) {
    UpdateMotorPosition();
  }
  // record previous position based on updated position
  UpdatePrevPosition();
  // check probability to bind/unbind
  if (UpdatePriors()) {
    // slow binding;
    // if binding or unbinding, don't diffuse
    if (rate_kinetics_ = SLOW) {
      return;
    }
    else if (rate_kinetics_ = FAST) {
      Diffuse()
    }
  }
  else if (bound_) {
    DiffuseAlongBond();
  }
  else {
    Diffuse();
  }
}

void Motor::Diffuse() {
  if (bound_) {
    // update motor step along bond
    DiffuseAlongBond();
  }
  else {
    // diffuse normally
  }
}

// update binding probabilities based on relative position to nearby bonds, then do rolls to determine if we bind or unbind
bool Motor::UpdatePriors() {
  // returns true if bound status is changed, else false

}
