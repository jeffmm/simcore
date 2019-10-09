#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink_species.hpp"
#include "params_parser.hpp"
#include "species_factory.hpp"
#include "output_manager.hpp"

class CrosslinkManager {
private:
  system_parameters *params_;
  OutputManager output_mgr_;
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
  void InitOutputs(bool reading_inputs = false, run_options *run_opts = nullptr);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void InitSpecies(sid_label &slab, ParamsParser &parser);
  void LoadCrosslinksFromCheckpoints(std::string run_name,
                                     std::string checkpoint_run_name);
  void ZeroDrTot();
  const double GetDrMax();
};

#endif
