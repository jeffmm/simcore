#include "crosslink_manager.h"

void CrosslinkManager::Init(system_parameters *params, 
    space_struct * space, MinimumDistance * mindist,
    std::vector<Object*> * objs) {
  params_ = params;
  space_ = space;
  mindist_ = mindist;
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
    update_ = true;
  }
  /* Check crosslink binding */
  if (gsl_rng_uniform_pos(rng_.r) <= xlink_concentration_*obj_volume_*k_on_*params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
    update_ = true;
  }
}

/* A crosslink binds to an object from solution */
void CrosslinkManager::BindCrosslink() {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  Crosslink xl;
  xlinks_singly_.push_back(xl);
  xlinks_singly_.back().Init(mindist_);
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
  //printf("%d \n", n_anchors_bound_);
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
  xlinks_singly_.back().SetSingly();
}

void CrosslinkManager::GetInteractors(std::vector<Object*> * ixors) {
  ixors->clear();
  std::vector<Object *> ix;
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++xlink) {
    ix.push_back(&(*xlink));
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void CrosslinkManager::UpdateCrosslinks() {
  /* First update bound crosslinks, then update binding */
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++xlink) {
    xlink->UpdateCrosslink();
  }
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end(); ++xlink) {
    xlink->UpdateCrosslink();
  }
  /* Check if we had any crosslinking events */
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end(); ++xlink) {
    if (xlink->IsDoubly()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_singly_.end()-1);
      xlinks_doubly_.push_back(*(xlinks_singly_.end()-1));
      xlinks_singly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_++;
    }
  }
  /* Check if we had any crosslinks break */
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end(); ++xlink) {
    if (!xlink->IsDoubly()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_doubly_.end()-1);
      xlinks_singly_.push_back(*(xlinks_doubly_.end()-1));
      xlinks_doubly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_--;
    }
  }

  /* Calculate binding of free crosslinks (in solution)
     and unbinding of bound crosslinks */
  CalculateBindingFree();
}

void CrosslinkManager::Clear() {
  xlinks_singly_.clear();
  xlinks_doubly_.clear();
}

void CrosslinkManager::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto it=xlinks_singly_.begin(); it!=xlinks_singly_.end(); ++it) {
    it->Draw(graph_array);
  }
  for (auto it=xlinks_doubly_.begin(); it!=xlinks_doubly_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToXlink(Object * xlink, Object * neighbor) {
  if (xlink->GetSID() != +species_id::crosslink) {
    error_exit("BindCrosslinkObj expected crosslink object, got generic object.");
  }
  Crosslink * xl = dynamic_cast<Crosslink*>(xlink);
  xl->AddNeighbor(neighbor);
}

