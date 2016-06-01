#ifndef _CYTOSCORE_FORCES_H_
#define _CYTOSCORE_FORCES_H_

#include "auxiliary.h"
#include "species.h"
#include "cell_list.h"
#include "minimum_distance.h"
#include "potential_manager.h"

class Forces {
  private:
    int n_dim_,
        n_periodic_;
    double dr_[3],
           contact1_[3],
           contact2_[3],
           dr_mag_,
           dr_mag2_,
           buffer_mag_,
           buffer_mag2_;
    space_struct *space_;
    std::vector<Simple*> simples_; 
    std::vector<cell_interaction> interactions_;
    CellList cell_list_;
    PotentialManager potentials_;
  public:
    void Init(space_struct *space, std::vector<SpeciesBase*> species, double cell_length);
    void UpdateCellList(std::vector<SpeciesBase*> species);
    void LoadSimples(std::vector<SpeciesBase*> species);
    void Interact();
    void CheckOverlap(std::vector<SpeciesBase*> species);
    void InitPotentials(std::vector<SpeciesBase*> species);
    void MinimumDistance(cell_interaction ix);
    interaction FirstInteraction(PotentialBase *pot);
    interaction SecondInteraction(PotentialBase *pot);
};

#endif // _CYTOSCORE_FORCES_H_
