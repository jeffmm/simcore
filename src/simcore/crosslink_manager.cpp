#include "crosslink_manager.h"

void CrosslinkManager::Init(system_parameters *params, std::vector<Object*> * objs) {
  params_ = params;
  objs_ = objs;
  /*TODO RNG should be initialized from a passed seed value */
  rng_.Init();
  k_on_ = params_->crosslink.k_on;
  k_off_ = params_->crosslink.k_off;
  xlink_concentration_ = params_->crosslink.concentration;
  obj_volume_ = 0;
  n_xlinks_ = 0;
  n_anchors_bound_ = 0;
  update_ = false;
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsVolume() {
  obj_volume_ = 0;
  for (auto it=objs_->begin(); it!=objs_->end(); ++it) {
    obj_volume_ += (*it)->GetVolume();
  }
  /* TODO, keep track of individual object binding probabilities,
   * which would end up being individual object volume/total
   * object volume. For now, will work correctly with equally-sized
   * objects. */
}

/* Whether to reinsert anchors into the interactors list */
bool CrosslinkManager::CheckUpdate() {
  if (update_) {
    update_ = false;
    return true;
  }
  return false;
}

void CrosslinkManager::CalculateBindingFree() {
  // Check crosslink unbinding
  if (gsl_rng_uniform_pos(rng_.r) <= k_off_*n_anchors_bound_*params_->delta) {
    /* Remove a random anchor from the system */
    UnbindCrosslink();
  }
  /* Check crosslink binding */
  if (gsl_rng_uniform_pos(rng_.r) <= xlink_concentration_*obj_volume_*k_on_*params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
  }
}

/* A crosslink binds to an object from solution */
void CrosslinkManager::BindCrosslink() {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  Crosslink xl;
  xlinks_singly_.push_back(xl);
  xlinks_singly_.back().Init(params_);
  /* Attach to random object in system */
  /* TODO Should weight probability of selecting object
     by object volume in the case of different sized
     objects */
  int i_obj = gsl_rng_uniform_int(rng_.r, objs_->size());
  xlinks_singly_.back().AttachObjRandom((*objs_)[i_obj]);
  /* Keep track of bound bound anchors, bound crosslinks, and
   * concentration of free crosslinks */
  n_anchors_bound_++;
  n_xlinks_++;
  xlink_concentration_ -= 1.0/space_->volume;
}

/* An unbinding event for a single, random anchor in the system */
void CrosslinkManager::UnbindCrosslink() {
  int i_anchor = gsl_rng_uniform_int(rng_.r, n_anchors_bound_);
  /* If i_anchor < the number of singly-bound anchors, unbind a
   * singly-bound anchor */
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
  if (xlinks_singly_.empty()) {
    return;
  }
  if (xlinks_singly_.size() > 1) {
    auto it = (xlinks_singly_.begin() + i_xlink);
    it->UnbindAnchor();
    std::iter_swap(xlinks_singly_.begin() + i_xlink, xlinks_singly_.end()-1);
    xlinks_singly_.pop_back();
  }
  else {
    xlinks_singly_.clear();
  }
  n_xlinks_--;
}

void CrosslinkManager::DoublyToSingly(int i_doubly) {
  if (xlinks_doubly_.size() > 1) {
    std::iter_swap(xlinks_doubly_.begin() + i_doubly, xlinks_doubly_.end()-1);
    xlinks_singly_.push_back(*(xlinks_doubly_.end()-1));
    xlinks_doubly_.pop_back();
  }
  else {
    xlinks_singly_.push_back(*(xlinks_doubly_.begin()));
    xlinks_doubly_.clear();
  }
}

void CrosslinkManager::GetInteractors(std::vector<Object*> * ixors) {
  ixors->clear();
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++xlink) {
    ixors->push_back(xlink->GetBoundPtr());
  }
}

void CrosslinkManager::UpdateCrosslinks() {
  /* First update bound crosslinks, then update binding */
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++xlink) {
    xlink->UpdateCrosslink();
  }
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end(); ++xlink) {
    xlink->UpdateCrosslink();
  }
  /* Calculate binding of free crosslinks (in solution)
     and unbinding of bound crosslinks */
  CalculateBindingFree();
}

void CrosslinkManager::Clear() {
  xlinks_singly_.clear();
  xlinks_doubly_.clear();
}
