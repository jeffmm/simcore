#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink.h"

class CrosslinkManager {
  private:
    bool update_;
    int n_xlinks_,
        n_anchors_bound_;
    MinimumDistance * mindist_;
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
    void CalculateBindingFree();
    void BindCrosslink();
    void UnbindCrosslink();
    void DoublyToSingly(int i_doubly);
    void RemoveCrosslink(int i_xlink);
  public:
    void Init(system_parameters *params, space_struct * space,
        MinimumDistance * mindist, std::vector<Object*> * objs);
    void GetInteractors(std::vector<Object*> * ixors);
    void UpdateCrosslinks();
    void UpdateObjsVolume();
    bool CheckUpdate();
    void Clear();
    void Draw(std::vector<graph_struct*> * graph_array);
    void BindCrosslinkObj(Object * obj);
    void AddNeighborToXlink(Object * xlink, Object * neighbor);
};

#endif
