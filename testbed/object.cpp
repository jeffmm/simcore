#include "object.h"
void generate_random_unit_vector(double * vec, unsigned int const n_dim, gsl_rng * r);
void get_random_coordinate(double * pos, unsigned int const n_dim, double radius, boundary_type btype, gsl_rng * r);


/****************************
** Object member functions **
****************************/
Object::Object() {
  oid_ = _next_oid_++; 
  n_dim_ = _n_dim_;
  std::fill(position_,position_+3,0.0);
  std::fill(orientation_,orientation_+3,0.0);
  length_ = 0;
  diameter_ = 0;
  rng_.Init(_seed_);
  _seed_ = gsl_rng_get(rng_.r);
}

unsigned int Object::_next_oid_ = 0;
unsigned int Object::_n_dim_ = 0;
long Object::_seed_ = 7777777;
unsigned int Object::GetOID() {
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
void Object::Report() {
  fprintf(stderr, "    OID: %d\n",GetOID());
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n",position_[0],position_[1],position_[2]);
  fprintf(stderr, "      u: {%2.2f %2.2f %2.2f}\n",orientation_[0],orientation_[1],orientation_[2]);
  fprintf(stderr, "      l: %2.2f\n",length_);
  fprintf(stderr, "      d: %2.2f\n",diameter_);
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
Bond * Site::GetOtherBond(unsigned int bond_oid) {
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
  sites_.push_back(s);
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
// Default d=1
void Mesh::AddRandomBondAnywhere(double l, double d) {
  if (n_sites_ < 1) {
    InitRandomSite(d);
  }
  int i_site=gsl_rng_uniform_int(rng_.r,n_sites_);
  AddRandomBondToSite(l, i_site);
}
void Mesh::AddRandomBondToSite(double l, int i_site) {
  if (i_site > n_sites_) {
    std::cerr << "Error!\n";
  }
  double const d = sites_[i_site].GetDiameter();
  double const * const pos0 = sites_[i_site].GetPosition();
  double pos[3];
  generate_random_unit_vector(pos, n_dim_, rng_.r);
  for (int i=0; i<n_dim_; ++i) {
    pos[i] *= l;
    pos[i] += pos0[i];
  }
  InitSiteAt(pos,d);

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
  ReportSites();
  ReportBonds();
}
// Default d=1
void Mesh::InitRandomSite(double d) {
  double pos[3];
  get_random_coordinate(pos,params_->n_dim,params_->system_radius,params_->boundary,rng_.r);
  InitSiteAt(pos,d);
}
