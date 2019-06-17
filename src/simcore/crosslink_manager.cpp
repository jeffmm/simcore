#include "crosslink_manager.h"

void CrosslinkManager::UpdateObjsVolume() {
  obj_volume_ = 0;
  for (auto it=objs_->begin(); it!=objs_->end(); ++it) {
    obj_volume_ += (*it_)->GetVolume()
  }
}

void CrosslinkManager::CalculateBinding() {
  if (++n_step_ % 100 != 0) return;
  double vol = GetVolume();
  // Check crosslink unbinding
  if (gsl_rng_uniform_pos(rng_.r) <= k_off_*n_xlinks_bound_*100*delta_) {
    // Remove a random anchor from the system
    UnbindCrosslink();
  }
  // Check crosslink binding
  if (gsl_rng_uniform_pos(rng_.r) <= xlink_concentration_*vol*k_on_*100*delta_) {
    // Create a new crosslink and bind an anchor to a random object in the system
    BindCrosslink();
  }
}

void CrosslinkManager::BindCrosslink() {
  Crosslink xl;
  xlinks_.push_back(xl);
  xlinks_.back().Init();
  int i_bond = gsl_rng_uniform_int(rng_.r,n_bonds_);
  xlinks_.back().SetMeshID(GetMeshID());
  xlinks_.back().AttachBondRandom(&bonds_[i_bond],i_bond*bond_length_);
  n_anchors_bound_++;
  xlink_concentration_ -= 1.0/space_->volume;
}

void CrosslinkManager::UnbindCrosslink() {
  int i_anchor = gsl_rng_uniform_int(rng_.r,n_anchors_bound_);
  int n_singly = xlinks_singly_.size();
  if (i_anchor < n_singly) {
    RemoveCrosslink(i_anchor);
  }
  else {
    i_anchor -= n_singly;
    int i_doubly = i_anchor/2;
    int i_which = i_anchor%2;
    xlinks_doubly_[i_doubly].UnbindAnchor(i_which);
    DoublyToSingly(i_doubly);
  }
  n_anchors_bound_--;
  xlink_concentration_ += 1.0/space_->volume;
}

void CrosslinkManager::RemoveCrosslink(int i_xlink) {
  if (xlinks_singly_.size() > 1) {
    auto it = (xlinks_singly_.begin() + i_xlink);
    it->UnbindAnchor();
    std::iter_swap(xlinks_singly_.begin() + i_xlink, xlinks_singly_.end()-1);
    xlinks_singly_.pop_back();
  }
  else {
    xlinks_singly_.clear();
  }
}

void CrosslinkManager::DoublyToSingly(int i_doubly) {
  if (xlinks_doubly_.size() > 1) {
    std::iter_swap(xlinks_doubly_.begin() + i_doubly, xlinks_doubly_.end()-1);
    xlinks_singly_.push_back(xlinks_doubly_.end()-1);
    xlinks_doubly_.pop_back();
  }
  else {
    xlinks_singly_.push_back(xlinks_doubly_.begin());
    xlinks_doubly_.clear();
  }
}

std::vector<Object*> CrosslinkManager::GetInteractors() {
  std::vector<Object*> interactors;
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++it) {
    interactors.push_back(xlink->GetBoundPtr());
  }
  return interactors;
}
