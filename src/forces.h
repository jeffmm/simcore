#ifndef _CYTOSCORE_FORCES_H_
#define _CYTOSCORE_FORCES_H_

#include "auxiliary.h"
#include "species.h"
#include "cell_list.h"

class Forces {
  private:
    int n_dim_,
        n_periodic_;
    space_struct *space_;
    std::vector<Simple*> simples_; 
    std::vector<cell_interaction> interactions_;
    CellList cell_list_;
  public:
    void Init(space_struct *space, double cell_length) {
      space_=space;
      n_dim_ = space->n_dim;
      n_periodic_ = space->n_periodic;
      cell_list_.Init(n_dim_, n_periodic_, cell_length, space->radius);
    }
    void UpdateCellList(std::vector<SpeciesBase*> species) {
      LoadSimples(species);
      interactions_.clear();
      cell_list_.LoadSimples(simples_);
      interactions_ = cell_list_.GetInteractions();
    }
    void LoadSimples(std::vector<SpeciesBase*> species) {
      simples_.clear();
      for (auto it=species.begin(); it!=species.end(); ++it) {
        std::vector<Simple*> sim_vec = (*it)->GetSimples();
        simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
      }
    }
    void CheckOverlap();
};

#endif // _CYTOSCORE_FORCES_H_
