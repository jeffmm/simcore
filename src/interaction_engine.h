#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_
 

//#include "cell_list.h"
#include "species.h"
#include "auxiliary.h"
#include "potential_manager.h"
#include "minimum_distance.h"
#include "particle_tracker.h"
#include "struct_analysis.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {

  private:
    double stress_[9],
           dr_update_;
    bool overlap_,
         no_interactions_,
         no_init_ = true,
         processing_ = false;
    int n_dim_,
        n_periodic_,
        i_update_,
        n_update_,
        n_objs_,
        n_thermo_,
        static_pnumber_;
    system_parameters *params_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;
    std::vector<ix_pair> nlist_;

    MinimumDistance mindist_;
    StructAnalysis struct_analysis_;

    std::vector<pair_interaction> pair_interactions_;
    std::vector<boundary_interaction> boundary_interactions_;
    std::vector<Object*> interactors_;
    //CellList clist_;
    ParticleTracker ptracker_;
    PotentialManager potentials_;

    void CheckUpdate();
    void UpdateInteractors();
    void UpdateInteractions();
    void ProcessPairInteraction(std::vector<pair_interaction>::iterator pix);
    void ProcessBoundaryInteraction(std::vector<boundary_interaction>::iterator bix);
    void CalculatePairInteractionsMP();
    void CalculateBoundaryInteractionsMP();
    void CalculatePairInteractions();
    void CalculateBoundaryInteractions();
    void ApplyPairInteractions();
    void ApplyBoundaryInteractions();
    double GetDrMax();
    void ZeroDrTot();
    int CountSpecies();

  public:
    InteractionEngine() {}
    void Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species, 
              space_struct *space,
              bool processing=false);
    void Interact();
    void CalculatePressure();
    bool CheckOverlap(std::vector<Object*> ixs);
    bool CheckBoundaryConditions(std::vector<Object*> ixs);
    void AddInteractors(std::vector<Object*> ixs);
    void Reset();
    void Clear();
    void StructureAnalysis();
    void CalculateStructure();
    void ForceUpdate();
};

#endif
