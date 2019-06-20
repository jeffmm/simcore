#include "crosslink.h"

void Crosslink::Init(system_parameters * params) {
  params_ = params;
}

void Crosslink::UnbindAnchor(bool second) {
  // If singly bound, we have to determine which anchor is bound
  if (!doubly_bound_) {
    if (anchors_.first.IsBound()) {
      anchors_.first.Clear();
    }
    else {
      anchors_.second.Clear();
    }
  }
  if (second) {
    anchors_.second.Clear();
  }
  else {
    anchors_.first.Clear();
  }
}

// Get pointer to anchor that is bound. Should only be used for singly bound xlinks.
Anchor * Crosslink::GetBoundPtr() {
  if (doubly_bound_) {
    warning("GetBoundPtr() called on doubly bound crosslink. This function"
        "should only be called on crosslinks with a free anchor!");
  }
  if (anchors_.first.IsBound()) {
    return &(anchors_.first);
  }
  else {
    return &(anchors_.second);
  }
}

void Crosslink::UpdateCrosslink() {
  if (doubly_bound_) {
    anchors_.first.UpdatePosition();
    anchors_.second.UpdatePosition();
  }
  else {
    if (anchors_.first.IsBound()) {
      anchors_.first.UpdatePosition();
    }
    else {
      anchors_.second.UpdatePosition();
    }
  }
  // TODO Apply tether forces, torques
}

/* Attach a crosslink anchor to object in a random fashion. Currently only bonds are used, but can be extended to sites (sphere-like objects) */
void Crosslink::AttachObjRandom(Object * obj) {
  /* Attaching to random bond implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj->GetType() == +obj_type::bond) {
    anchors_.first.AttachObjRandom(obj);
  }
  else {
    /* TODO: add binding to sphere or site-like objects */
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
}

