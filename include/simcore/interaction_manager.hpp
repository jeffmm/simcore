#ifndef _SIMCORE_INTERACTION_MANAGER_H_
#define _SIMCORE_INTERACTION_MANAGER_H_

//#include "cell_list.hpp"
#include "auxiliary.hpp"
#include "cell_list.hpp"
#include "crosslink_manager.hpp"
#include "minimum_distance.hpp"
#include "potential_manager.hpp"
#include "species.hpp"
#include "struct_analysis.hpp"

typedef std::vector<Interaction>::iterator ix_iterator;

class InteractionManager {
 private:
  double stress_[9];
  double dr_update_;
  bool overlap_;
  bool no_interactions_;
  bool no_boundaries_;
  bool no_init_ = true;
  bool processing_ = false;
  bool in_out_flag_ = false;
  int n_dim_;
  int n_periodic_;
  int n_objs_;
  int n_thermo_;
  int static_pnumber_;
  int *i_step_;
  int n_interactions_;
  int i_update_ = 0;
  system_parameters *params_;
  space_struct *space_;
  std::vector<SpeciesBase *> *species_;

  MinimumDistance mindist_;
  StructAnalysis struct_analysis_;

  std::vector<Interaction> pair_interactions_;
  std::vector<Interaction> boundary_interactions_;
  std::vector<Object *> ix_objects_;
  std::vector<Object *> interactors_;
  CellList clist_;
  PotentialManager potentials_;
  CrosslinkManager xlink_;

  const bool CheckSpeciesInteractorUpdate() const;
  void CheckUpdateXlinks();
  void CheckUpdateInteractions();
  void UpdateInteractors();
  void UpdateInteractions();
  void UpdatePairInteractions();
  void UpdateBoundaryInteractions();
  void ProcessPairInteraction(ix_iterator ix);
  void ProcessBoundaryInteraction(ix_iterator ix);
  void CalculatePairInteractions();
  void CalculateBoundaryInteractions();
  void ApplyPairInteractions();
  void ApplyBoundaryInteractions();
  double GetDrMax();
  void ZeroDrTot();
  int CountSpecies();
  bool CheckBondAnchorPair(Object *anchor, Object *bond);

 public:
  InteractionManager() {}
  void Init(system_parameters *params, std::vector<SpeciesBase *> *species,
            space_struct *space, bool processing = false);
  void InitInteractions();
  void Interact();
  void CalculatePressure();
  bool CheckOverlap(std::vector<Object *> &ixs);
  bool CheckBoundaryConditions(std::vector<Object *> &ixs);
  void AddInteractors(std::vector<Object *> &ixs);
  void Reset();
  void Clear();
  void StructureAnalysis();
  void CalculateStructure();
  void ForceUpdate();
  void CheckUpdateObjects();
  void DrawInteractions(std::vector<graph_struct *> &graph_array);
  void WriteOutputs();
  void InitOutputs(bool reading_inputs = false,
                   run_options *run_opts = nullptr);
  void ReadInputs();
  void ResetCellList();
  void InitCrosslinkSpecies(sid_label &slab, ParamsParser &parser,
                            unsigned long seed);
  void LoadCrosslinksFromCheckpoints(std::string run_name,
                                     std::string checkpoint_run_name);
  void InsertCrosslinks();
  void BeginWithCrosslinks();
  void BindCrosslink();
  void PairBondCrosslinks();
};

#endif
