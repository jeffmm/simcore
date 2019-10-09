#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "params_parser.hpp"
#include "crosslink_species.hpp"
#include "species_factory.hpp"

class CrosslinkManager {
private:
  system_parameters *params_;
  double obj_volume_;
  bool update_;
  std::vector<CrosslinkSpecies *> xlink_species_;
  std::vector<Object *> *objs_;
  space_struct *space_;

public:
  void Init(system_parameters *params, space_struct *space,
            std::vector<Object *> *objs);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdateCrosslinks();
  void UpdateObjsVolume();
  bool CheckUpdate();
  void Clear();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void WriteOutputs();
  void ReadInputs();
  void InitOutputs(bool reading_inputs = false, bool reduce_flag = false,
                   bool with_reloads = false);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void InitSpecies(sid_label &slab, ParamsParser &parser);
  void ZeroDrTot();
  const double GetDrMax();
};

#endif
