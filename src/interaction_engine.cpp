
#include "interaction_engine.h"

void InteractionEngine::Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species, 
              space_struct *space) {
  // Set up pointer structures
  params_ = params;
  species_ = species;
  space_ = space;

  // Initialize owned structures
  n_dim_ = params_->n_dim;
  n_periodic_ = params_->n_periodic;
  n_update_ = params_->n_update_cells;
  i_update_ = 0;
  clist_.Init(n_dim_,n_periodic_,params_->cell_length,space_->radius);
  wca_.Init(params_);
  UpdateSimples();
  UpdateInteractions();
}

void InteractionEngine::UpdateSimples() {
  simples_.clear();
  for (auto spec_it= species_->begin(); spec_it!=species_->end(); ++spec_it) {
    std::vector<Simple*> simps = (*spec_it)->GetSimples();
    simples_.insert(simples_.end(),simps.begin(),simps.end());
  }
}

void InteractionEngine::UpdateInteractions() {
  clist_.LoadSimples(simples_);
  pair_interactions_ = clist_.GetPairInteractions();
}

void InteractionEngine::CheckUpdate() {
  if ((++i_update_) % n_update_ == 0)
    UpdateInteractions();
}


/****************************************
  INTERACT: Loop through interactions and
    apply WCA potential (for now)
*****************************************/

void InteractionEngine::Interact() {
  // Check if we need to update cell list
  CheckUpdate();
  // Loop through interactions
  for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix)
    ProcessInteraction(pix);
}

void InteractionEngine::ProcessInteraction(std::vector<pair_interaction>::iterator pix) {
  // Avoid certain types of interactions
  auto obj1 = pix->first.first;
  auto obj2 = pix->first.second;
  Interaction *ix = &(pix->second);
  if (obj1->GetRID() == obj2->GetRID())  return;
  if (obj1->GetCID() == obj2->GetCID())  return;
  if (obj1->GetOID() == obj2->GetOID()) 
    error_exit("ERROR! Object %d attempted self-interaction!\n",
        obj1->GetOID());

  //interactionmindist imd;
  MinimumDistance(obj1,obj2,ix,n_dim_,n_periodic_,space_);

  // XXX Only calculate WCA potential for now
  if (ix->dr_mag2 > wca_.GetRCut2())  return;
  double force[3];
  double torque[3];
  double pot_en;
  wca_.CalcPotential(ix);
  // Apply forces
  obj1->AddForce(ix->force);
  obj2->SubForce(ix->force);
  // XXX Assume object is extended, add torques
  cross_product(ix->contact1,ix->force,ix->t1,3);
  obj1->AddTorque(ix->t1);
  cross_product(ix->contact2,ix->force,ix->t2,3);
  obj2->SubTorque(ix->t2);
  // Add potential
  obj1->AddPotential(ix->pote);
  obj2->AddPotential(ix->pote);
}


