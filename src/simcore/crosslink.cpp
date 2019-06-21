#include "crosslink.h"

Crosslink::Crosslink() : Object() {
  SetSID(species_id::crosslink);
}

void Crosslink::Init(MinimumDistance * mindist) {
  mindist_ = mindist;
  length_ = -1;
  diameter_ = params_->crosslink.diameter;
  color_ = params_->crosslink.color;
  draw_ = draw_type::_from_string(params_->crosslink.draw_type.c_str());
  doubly_bound_ = false;
  /* TODO generalize crosslinks to more than two anchors */
  Anchor anchor1, anchor2;
  anchors_.push_back(anchor1);
  anchors_.push_back(anchor2);
  anchors_[0].Init();
  anchors_[1].Init();
}

void Crosslink::UpdatePosition() {
  SetPosition(anchors_[0].GetPosition());
  SetScaledPosition(anchors_[0].GetScaledPosition());
  SetOrientation(anchors_[0].GetOrientation());
  length_ = 0;
}

void Crosslink::UnbindAnchor(bool second) {
  /* If singly bound, just clear both */
  if (!doubly_bound_) {
    anchors_[0].Clear();
    anchors_[1].Clear();
    return;
  }
  if (second) {
    anchors_[1].Clear();
  }
  /* When singly bound, only consider first anchor populated */
  else {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
  }
}

/* Get pointer to anchor that is bound. Should only be used for singly bound
 * xlinks. */
Anchor * Crosslink::GetBoundPtr() {
  if (doubly_bound_) {
    warning("GetBoundPtr() called on doubly bound crosslink. This function"
        "should only be called on crosslinks with a free anchor!");
  }
  /* Only first anchor should be bound in the case of singly bound */
  return &(anchors_[0]);
}

void Crosslink::UpdateCrosslink() {
  anchors_[0].UpdatePosition();
  anchors_[1].UpdatePosition();
  if (doubly_bound_) {
    CalculateTetherForces();
  }
  else {
    UpdatePosition();
  }
}

void Crosslink::CalculateTetherForces() {
  mindist_->ObjectObject(&(anchors_[0]), &(anchors_[1]), 
      &ix);
  /* Check stretch of tether. No penalty for having a stretch < rest_length. ie
   * the spring does not resist compression. */
  length_ = sqrt(ix.dr_mag2);
  double stretch = length_ - rest_length_;
  std::fill(force_, force_+3, 0.0);
  for (int i=0; i<params_->n_dim; ++i) {
    orientation_[i] = ix.dr[i]/length_;
    position_[i] = ix.midpoint[i];
  }
  if (stretch > 0) {
    for (int i=0; i<params_->n_dim; ++i) {
      force_[i] = k_spring_ * stretch * orientation_[i];
    }
    anchors_[0].AddForce(force_);
    anchors_[1].SubForce(force_);
  }
  /* TODO Apply torques on crosslinks if necessary */
  
}

/* Attach a crosslink anchor to object in a random fashion. Currently only
 * bonds are used, but can be extended to sites (sphere-like objects) */
void Crosslink::AttachObjRandom(Object * obj) {
  /* Attaching to random bond implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj->GetType() == +obj_type::bond) {
    anchors_[0].AttachObjRandom(obj);
  }
  else {
    /* TODO: add binding to sphere or site-like objects */
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
}

void Crosslink::Draw(std::vector<graph_struct*> * graph_array) {
  /* Draw anchors */
  anchors_[0].Draw(graph_array);
  anchors_[1].Draw(graph_array);
  /* Draw tether */
  if (doubly_bound_ && length_ > 0) {
    std::copy(position_, position_+3, g_.r);
    std::copy(orientation_, orientation_+3, g_.u);
    g_.color = color_;
    if (params_->graph_diameter > 0) {
      g_.diameter = params_->graph_diameter;
    }
    else {
      g_.diameter = diameter_;
    }
    g_.length = length_;
    g_.draw = draw_;
    graph_array->push_back(&g_);
  }
}

void Crosslink::SetDoubly() {
  doubly_bound_ = true;
}

void Crosslink::SetSingly() {
  doubly_bound_ = false;
}
