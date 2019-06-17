#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.h"
#include "anchor.h"
#include "minimum_distance.h"

class Crosslink {
  private:
    bool doubly_bound_;

    double k_on_,
           k_off_,
           k_spring_,
           rest_length_;
    std::pair<Anchor, Anchor> anchors_;
  public:
    void UnbindAnchor(bool second);
    Anchor * GetBoundPtr();
};

#endif
