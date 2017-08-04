#include "bond.h"

/**************************
** Site member functions **
**************************/
void Site::AddBond(Bond * bond, directed_type dir) {
  bonds_.push_back(std::make_pair(bond,dir));
  n_bonds_++;
}
Bond * Site::GetBond(int i) {
  if (i<0 || i >= bonds_.size()) {
    std::cerr << "ERROR! Requested adjacent bond out of bounds!\n";
  }
  return bonds_[i].first;
}
Bond * Site::GetOtherBond(int bond_oid) {
  return GetOtherDirectedBond(bond_oid).first;
}

directed_bond Site::GetOtherDirectedBond(int bond_oid) {
  if (n_bonds_ == 1) {
    // there are no other bonds
    return std::make_pair(nullptr,NONE);
  }
  int i_bond = gsl_rng_uniform_int(rng_.r, n_bonds_-1);
  if (bonds_[i_bond].first->GetOID() != bond_oid) {
    return bonds_[i_bond];
  }
  return bonds_[n_bonds_-1];
}

void Site::Report() {
  fprintf(stderr, "  Site:\n");
  Object::Report();
}
void Site::ReportBonds() {
  Report();
  fprintf(stderr,"    Reporting bonds:\n");
  for (std::vector<directed_bond>::iterator it=bonds_.begin(); it!=bonds_.end(); ++it) {
    (*it).first->Report();
  }
}


