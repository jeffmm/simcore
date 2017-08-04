#include "mesh.h"

/**************************
** Mesh member functions **
**************************/
Mesh::Mesh(long seed) : Object(seed) {
  n_sites_ = n_bonds_ = 0;
}

void Mesh::Reserve(int n_bonds) {
  sites_.reserve(n_bonds+1);
  bonds_.reserve(n_bonds);
  n_bonds_max_ = n_bonds;
}

Site * Mesh::GetSite(int i) {
  return &(sites_[i]);
}
Bond * Mesh::GetBond(int i) {
  return &(bonds_[i]);
}
void Mesh::AddSite(Site s) {
  if (n_sites_ == n_bonds_max_+1) {
    std::cerr << "ERROR! Attempting to add site beyond allocated maximum.\n";
    exit(1);
  }
  sites_.push_back(s);
  n_sites_++;
}
void Mesh::AddBond(Bond b) {
  if (n_bonds_ == n_bonds_max_) {
    std::cerr << "ERROR! Attempting to add bond beyond allocated maximum.\n";
    exit(1);
  }
  bonds_.push_back(b);
  n_bonds_++;
}

void Mesh::InitSiteAt(double * pos, double d) {
  Site s(gsl_rng_get(rng_.r));
  s.SetPosition(pos);
  s.SetDiameter(d);
  AddSite(s);
}
void Mesh::InitBondAt(double * pos, double * u, double l, double d) {
  Site s1(gsl_rng_get(rng_.r));
  Site s2(gsl_rng_get(rng_.r));
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
  Bond b(gsl_rng_get(rng_.r));
  AddBond(b);
  bonds_[n_bonds_-1].Init(&sites_[n_sites_-2],&sites_[n_sites_-1]);
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
  double pos[3] = {0,0,0};
  generate_random_unit_vector(pos, n_dim_, rng_.r);
  for (int i=0; i<n_dim_; ++i) {
    pos[i] = pos0[i] + l * pos[i];
  }
  InitSiteAt(pos,d);
  Bond b(gsl_rng_get(rng_.r));
  AddBond(b);
  bonds_[n_bonds_-1].Init(&sites_[i_site],&sites_[n_sites_-1]);
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
  Bond b(gsl_rng_get(rng_.r));
  AddBond(b);
  bonds_[n_bonds_-1].Init(&sites_[i_site],&sites_[n_sites_-1]);
}
void Mesh::ReportSites() {
  for (std::vector<Site>::iterator it=sites_.begin(); it!=sites_.end(); ++it) {
    it->Report();
    std::cerr << &(*it) << "\n";
  }
}
void Mesh::ReportBonds() {
  for (std::vector<Bond>::iterator it=bonds_.begin(); it!=bonds_.end(); ++it) {
    it->Report();
    std::cerr << &(*it) << "\n";
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
void Mesh::SubReport() {
  for (std::vector<Bond>::iterator it=bonds_.begin(); it!=bonds_.end(); ++it) {
    it->ReportSites();
  }
  for (std::vector<Site>::iterator it=sites_.begin(); it!=sites_.end(); ++it) {
    it->ReportBonds();
  }
}
void Mesh::SetBondLength(double l) {
  bond_length_ = l;
}
void Mesh::Draw(std::vector<graph_struct*> * graph_array) {
  for (std::vector<Bond>::iterator it = bonds_.begin(); it!=bonds_.end(); ++it) {
    it->Draw(graph_array);
  }
}

