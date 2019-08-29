#ifndef _SIMCORE_CROSSLINK_MANAGER_H_
#define _SIMCORE_CROSSLINK_MANAGER_H_

#include "crosslink.hpp"

class CrosslinkManager {
 private:
  bool update_;
  int n_xlinks_, n_anchors_bound_, n_spec_, n_checkpoint_, spec_flag_,
      checkpoint_flag_;
  MinimumDistance *mindist_;
  std::string checkpoint_file_;
  double obj_volume_, xlink_concentration_, k_on_, k_off_;
  RNG rng_;
  space_struct *space_;
  LookupTable lut_;
  std::vector<Crosslink> xlinks_doubly_;
  std::vector<Crosslink> xlinks_singly_;
  std::vector<Object *> *objs_;
  std::fstream ispec_file_;
  std::fstream ospec_file_;
  system_parameters *params_;
  void CalculateBindingFree();
  void BindCrosslink();
  void UnbindCrosslink();
  void DoublyToSingly(int i_doubly);
  void RemoveCrosslink(int i_xlink);

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
  void AddNeighborToXlink(Object *xlink, Object *neighbor);
  void WriteOutputs();
  void ReadInputs();
  void InitOutputs(bool reading_inputs = false, bool reduce_flag = false,
                   bool with_reloads = false);
};

#endif
