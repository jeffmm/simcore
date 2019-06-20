#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.h"
#include "anchor.h"
#include "minimum_distance.h"

class Crosslink {
  private:
    system_parameters * params_;
    int mesh_id_;
    bool doubly_bound_;
    double k_on_,
           k_off_,
           k_spring_,
           k_align_,
           rest_length_;
    std::pair<Anchor, Anchor> anchors_;

  public:
    Crosslink() {}
    void Init(system_parameters *params);
    int const GetMeshID() const;
    void SetMeshID(int mid);
    void UpdateCrosslink();
    void UnbindAnchor(bool second=false);
    void AttachObjRandom(Object * obj);
    Anchor * GetBoundPtr();
};

#endif
