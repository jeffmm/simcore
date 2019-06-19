#include "crosslink.h"

void Crosslink::UnbindAnchor(bool second=false) {
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
