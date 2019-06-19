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
           k_align_,
           rest_length_;
    std::pair<Anchor, Anchor> anchors_;

  public:
    void Init(system_parameters *params, long seed);
    void UpdateCrosslink();
    void UnbindAnchor(bool second=false);
    Anchor * GetBoundPtr();
};

#endif
