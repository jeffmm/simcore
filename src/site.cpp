#include "bond.h"

/**************************
** Site member functions **
**************************/
Site::Site() : Object() {
  n_bonds_ = 0;
  std::fill(tangent_,tangent_+3,0.0);
  std::fill(random_force_,random_force_+3,0.0);
}
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

directed_bond Site::GetDirectedBond(int i) {
  if (i<0 || i >= bonds_.size()) {
    std::cerr << "ERROR! Requested adjacent bond out of bounds!\n";
  }
  return bonds_[i];
}

directed_bond Site::GetOutgoingBond() {
  for (auto it=bonds_.begin(); it!=bonds_.end(); ++it) {
    if ((*it).second == OUTGOING ) {
      return *it;
    }
  }
  return std::make_pair(nullptr,NONE);
}

void Site::RemoveBond(int bond_oid) {
  for (int i_bond=0; i_bond<n_bonds_;++i_bond) {
    if (bonds_[i_bond].first->GetOID() == bond_oid) {
      bonds_.erase(bonds_.begin()+i_bond);
      n_bonds_--;
      return;
    }
  }
}

void Site::RemoveOutgoingBonds() {
  for (int i_bond=0; i_bond<n_bonds_;) {
    if (bonds_[i_bond].second == OUTGOING) {
      bonds_.erase(bonds_.begin()+i_bond);
      n_bonds_--;
    }
    else {
      ++i_bond;
    }
  }
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

void Site::CalcTangent() {
  if (n_bonds_ == 1) {
    double const * const u = bonds_[0].first->GetOrientation();
    std::copy(u,u+3,tangent_);
  }
  else if (n_bonds_ == 2) {
    double const * const u1 = bonds_[0].first->GetOrientation();
    double const * const u2 = bonds_[1].first->GetOrientation();
    for (int i=0;i<n_dim_;++i) {
      tangent_[i] = u1[i] + u2[i];
    }
    normalize_vector(tangent_,n_dim_);
  }
  else {
    std::fill(tangent_,tangent_+3,0.0);
  }
}
void Site::SetRandomForce(double * f_rand) {
  std::copy(f_rand,f_rand+3,random_force_);
}
double const * const Site::GetTangent() {
  return tangent_;
}
double const * const Site::GetRandomForce() {
  return random_force_;
}
void Site::AddRandomForce() {
  for (int i=0; i<n_dim_; ++i) {
    force_[i] += random_force_[i];
  }
}
bool Site::HasNeighbor(int other_oid) {
  if (n_bonds_ == 0) {
    return false;
  }
  for (db_iterator db = bonds_.begin(); db!=bonds_.end(); ++db) {
    if (db->first->GetOID() == other_oid) {
      return true;
    }
  }
  return false;
}

//void Site::Draw(std::vector<graph_struct*> * graph_array) {
  //std::copy(scaled_position_,scaled_position_+3, g_.r);
  //for (int i=space_->n_periodic; i<n_dim_; ++i) {
    //g_.r[i] = position_[i];
  //}
  //std::copy(orientation_, orientation_+3, g_.u);
  //g_.color = color_;
  //g_.diameter = diameter_;
  //g_.length = length_;
  //g_.draw = draw_;
  //graph_array->push_back(&g_);
//}

