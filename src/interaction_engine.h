#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_
 

#include "cell_list.h"
#include "species.h"
#include "auxiliary.h"

#include "wca.h"
#include "interaction.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {

  private:
    double virial_;
    int n_dim_,
        n_periodic_,
        i_update_,
        n_update_,
        n_objs_,
        virial_time_avg_;
    system_parameters *params_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;

    std::vector<pair_interaction> pair_interactions_;
    std::vector<Simple*> simples_;
    CellList clist_;
    WCA wca_;

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
};

#endif
