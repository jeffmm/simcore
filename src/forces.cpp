#include "forces.h"

void Forces::Init(space_struct *space, std::vector<SpeciesBase*> species, double cell_length) {
  space_=space;
  n_dim_ = space->n_dim;
  n_periodic_ = space->n_periodic;
  cell_list_.Init(n_dim_, n_periodic_, cell_length, space->radius);
  CheckOverlap(species);
}
void Forces::UpdateCellList(std::vector<SpeciesBase*> species) {
  LoadSimples(species);
  interactions_.clear();
  cell_list_.LoadSimples(simples_);
  interactions_ = cell_list_.GetInteractions();
}
void Forces::LoadSimples(std::vector<SpeciesBase*> species) {
  simples_.clear();
  for (auto it=species.begin(); it!=species.end(); ++it) {
    std::vector<Simple*> sim_vec = (*it)->GetSimples();
    simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
  }

}
void Forces::Interact() {
  for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
    if (it->first->GetOID() == it->second->GetOID()) {
      error_exit("ERROR: Got self-interaction.\n");
    }
    if (it->first->GetCID() == it->second->GetCID()) continue;

  }
}
void Forces::CheckOverlap(std::vector<SpeciesBase*> species) {
  bool overlap = true;
  int num = 0;
  do {
    num++;
    overlap=false;
    UpdateCellList(species);
    for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
      if (it->first->GetCID() == it->second->GetCID()) continue;
      min_distance_point_point (n_dim_, n_periodic_, space_->unit_cell, 
          it->first->GetPosition(), it->first->GetScaledPosition(),
          it->second->GetPosition(), it->second->GetScaledPosition(),
          dr_, &dr_mag2_);
      double buffer_mag2 = 0.5*(it->first->GetDiameter() + 
                                it->second->GetDiameter());
      buffer_mag2=SQR(buffer_mag2);
      if (dr_mag2_ < buffer_mag2) {
        overlap = true;
        unsigned int const sid = it->second->GetSID();
        unsigned int const cid = it->second->GetCID();
        for (auto spec=species.begin(); spec!=species.end(); ++spec) {
          if ((*spec)->GetSID() == sid) {
            (*spec)->ReInit(cid);
          }
        }
        break;
      }
    }
    if (num > 10000)
      error_exit("ERROR: Too many overlaps detected. Check packing ratio for objects.\n");
  } while (overlap);
}



