#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink.hpp"

class CrosslinkManager {
private:
  bool update_;
  int n_xlinks_;
  int n_anchors_bound_;
  int n_spec_;
  int n_checkpoint_;
  int spec_flag_;
  int checkpoint_flag_;
  MinimumDistance *mindist_;
  std::string checkpoint_file_;
  double obj_volume_;
  double xlink_concentration_;
  double k_on_;
  double k_off_;
  double k_on_d_;
  double k_off_d_;
  RNG rng_;
  space_struct *space_;
  LookupTable lut_;
  std::vector<Crosslink> xlinks_;
  std::vector<Object *> *objs_;
  std::fstream ispec_file_;
  std::fstream ospec_file_;
  system_parameters *params_;
  void CalculateBindingFree();
  void BindCrosslink();
  void UpdateBoundCrosslinks();
  void ApplyCrosslinkTetherForces();
  Object *GetRandomObject();

  /* IO Functions */
  void WriteSpecs();
  void ReadSpecs();
  void WriteCheckpoints();
  void ReadCheckpoints();
  void InitOutputFiles();
  void InitSpecFile();
  void InitCheckpoints();
  bool InitSpecFileInputFromFile(std::string spec_file_name);
  void InitSpecFileInput();
  void LoadFromCheckpoints();

public:
  void Init(system_parameters *params, space_struct *space,
            MinimumDistance *mindist, std::vector<Object *> *objs);
  void GetInteractors(std::vector<Object *> *ixors);
  void UpdateCrosslinks();
  void UpdateObjsVolume();
  bool CheckUpdate();
  void Clear();
  void Draw(std::vector<graph_struct *> *graph_array);
  void BindCrosslinkObj(Object *obj);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void WriteOutputs();
  void ReadInputs();
  void InitOutputs(bool reading_inputs = false, bool reduce_flag = false,
                   bool with_reloads = false);
  void GetAnchorInteractors(std::vector<Object *> *ixors);
};

#endif
