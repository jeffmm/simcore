#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink_species.hpp"

class CrosslinkManager {
private:
  double obj_volume_;
  bool update_;
  std::vector<CrosslinkSpecies *> xlink_species_;
  std::vector<Object *> *objs_;
  void InitSpecies(system_parameters *params, space_struct *space,
                   MinimumDistance *mindist, std::vector<Object *> *objs);

public:
  // iENGINE USES
  void Init(system_parameters *params, space_struct *space,
            MinimumDistance *mindist, std::vector<Object *> *objs);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdateCrosslinks();
  void UpdateObjsVolume();
  bool CheckUpdate();
  void Clear();
  void Draw(std::vector<graph_struct *> *graph_array);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void WriteOutputs();
  void ReadInputs();
  void InitOutputs(bool reading_inputs = false, bool reduce_flag = false,
                   bool with_reloads = false);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
};

#endif
