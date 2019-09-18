#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_

//#include "cell_list.hpp"
#include "auxiliary.hpp"
#include "crosslink_manager.hpp"
#include "minimum_distance.hpp"
#include "particle_tracker.hpp"
#include "potential_manager.hpp"
#include "species.hpp"
#include "struct_analysis.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {
private:
  double stress_[9];
  double dr_update_;
  bool overlap_;
  bool no_interactions_;
  bool no_init_ = true;
  bool processing_ = false;
  bool in_out_flag_ = false;
  int n_dim_;
  int n_periodic_;
  int i_update_;
  int n_update_;
  int n_objs_;
  int n_thermo_;
  int static_pnumber_;
  int *i_step_;
  int n_interactions_;
  system_parameters *params_;
  space_struct *space_;
  std::vector<SpeciesBase *> *species_;
  std::vector<ix_pair> nlist_;

  MinimumDistance mindist_;
  StructAnalysis struct_analysis_;

  std::vector<pair_interaction> pair_interactions_;
  std::vector<boundary_interaction> boundary_interactions_;
  std::vector<Object *> ix_objects_;
  std::vector<Object *> interactors_;
  // CellList clist_;
  ParticleTracker ptracker_;
  PotentialManager potentials_;
  CrosslinkManager xlink_;

  void CheckUpdate();
  void UpdateInteractors();
  void UpdateInteractions();
  void ProcessPairInteraction(std::vector<pair_interaction>::iterator pix);
  void
  ProcessBoundaryInteraction(std::vector<boundary_interaction>::iterator bix);
  void CalculatePairInteractions();
  void CalculateBoundaryInteractions();
  void ApplyPairInteractions();
  void ApplyBoundaryInteractions();
  double GetDrMax();
  void ZeroDrTot();
  int CountSpecies();
  void PairBondCrosslinks();
  bool CheckBondAnchorPair(Object * anchor, Object * bond);

public:
  InteractionEngine() {}
  void Init(system_parameters *params, std::vector<SpeciesBase *> *species,
            space_struct *space, int *i_step, bool processing = false);
  void Interact();
  void CalculatePressure();
  bool CheckOverlap(std::vector<Object *> ixs);
  bool CheckBoundaryConditions(std::vector<Object *> ixs);
  void AddInteractors(std::vector<Object *> ixs);
  void Reset();
  void Clear();
  void StructureAnalysis();
  void CalculateStructure();
  void ForceUpdate();
  bool CountAndUpdate();
  void DrawInteractions(std::vector<graph_struct *> *graph_array);
  void WriteOutputs();
  void InitOutputs(bool reading_inputs = false, bool reduce_flag = false,
                   bool with_reloads = false);
  void ReadInputs();
};

#endif
