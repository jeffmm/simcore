#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink.h"

class CrosslinkManager {
  private:
    bool update_;
    int n_xlinks_,
        n_anchors_bound_;
    double obj_volume_,
           xlink_concentration_,
           k_on_,
           k_off_;
    RNG rng_;
    space_struct * space_;
    std::vector<Crosslink> xlinks_doubly_;
    std::vector<Crosslink> xlinks_singly_;
    std::vector<Object*> * objs_;
    system_parameters * params_;
    void UpdateObjsVolume();
    void CalculateBindingFree();
    void BindCrosslink();
    void UnbindCrosslink();
    void DoublyToSingly(int i_doubly);
    void RemoveCrosslink(int i_xlink);
  public:
    void Init(system_parameters *params, std::vector<Object*> * objs);
    void GetInteractors(std::vector<Object*> * ixors);
    void UpdateCrosslinks();
    bool CheckUpdate();
    void Clear();
};

#endif
