#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_
 

#include "cell_list.h"
#include "species.h"
#include "auxiliary.h"
#include "potential_manager.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {

  private:
    double stress_[9];
    bool overlap_,
         no_interactions_;
    int n_dim_,
        n_periodic_,
        i_update_,
        n_update_,
        n_objs_,
        n_thermo_;
    system_parameters *params_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;

    std::vector<pair_interaction> pair_interactions_;
    std::vector<Simple*> simples_;
    CellList clist_;
    PotentialManager potentials_;

    void CheckUpdate();
    void UpdateSimples();
    void UpdateInteractions();
    void ProcessInteraction(std::vector<pair_interaction>::iterator pix);
    void CalculateInteractionsMP();
    void CalculateInteractions();
    void ApplyInteractions();
    int CountSpecies();

  public:
    InteractionEngine() {}
    void Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species, 
              space_struct *space);
    void Interact();
    void CalculatePressure();
    bool CheckOverlap();
};

#endif
