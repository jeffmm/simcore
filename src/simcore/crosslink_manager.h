#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "species.h"
#include "auxiliary.h"
#include "minimum_distance.h"

class Crosslink {
  private:
    double k_on_,
           k_off_,
           bond_length_;
    std::vector<Anchor*> anchors_;
  public:
    void UnbindAnchor(int i_anchor) {
      auto it = anchors_.begin() + i_anchor;
      *(it)->Unbind();
      anchors_.erase(it);
    }
};

class CrosslinkManager {
  private:
    int n_xlinks_,
        n_anchors_bound_;
    double obj_volume_,
           xlink_concentration_,
           k_on_,
           k_off_;
    std::vector<Crosslink> xlinks_doubly_;
    std::vector<Crosslink> xlinks_singly_;
    std::vector<Object*> * objs_;
    system_parameters * params_;
    void UpdateObjsVolume();
    void CalculateBindingProbability();
    void BindCrosslink();
    void UnbindCrosslink();
    void DoublyToSingly();
    void RemoveCrosslink(int i_xlink);
};



#endif
