#ifndef _SIMCORE_FORCES_H_
#define _SIMCORE_FORCES_H_

#include "auxiliary.h"
#include "species.h"
#include "cell_list.h"
#include "minimum_distance.h"
#include "potential_manager.h"

// CJE stuff here
#include "force_brute.h"
#include "force_microcell.h"

class Forces {
  private:
    int n_dim_,
        n_periodic_;
    int draw_flag_;
    bool draw_;
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

    //ForceBase* test_force_;
    ForceMicrocell* test_microcell_f_;
  public:
    Forces() {
        //test_force_ = forceFactory<ForceBrute>();
        test_microcell_f_ = forceFactory<ForceMicrocell>();
    }
    ~Forces() {
        //delete(test_force_);
        delete(test_microcell_f_);
    }

    std::vector<graph_struct> draw_array_;
  public:
    void Init(space_struct *space, std::vector<SpeciesBase*> species, double cell_length, int draw_flag);
    void UpdateScheme(std::vector<SpeciesBase*> species);
    void UpdateCellList(std::vector<SpeciesBase*> species);
    void LoadSimples(std::vector<SpeciesBase*> species);
    void Interact();
    void InteractMP();
    void DumpAll();
    void CheckOverlap(std::vector<SpeciesBase*> species);
    void InitPotentials(std::vector<SpeciesBase*> species);
    void MinimumDistance(cell_interaction ix);
    void Draw(std::vector<graph_struct*> * graph_array);
    interaction FirstInteraction(PotentialBase *pot);
    interaction SecondInteraction(PotentialBase *pot);
};

#endif // _SIMCORE_FORCES_H_
