
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
  virial_time_avg_ = params_->virial_time_avg;
  i_update_ = 0;
  n_objs_ = CountSpecies();
  clist_.Init(n_dim_,n_periodic_,params_->cell_length,space_->radius);
  potentials_.InitPotentials(params_);
}

/****************************************
  INTERACT: Loop through interactions and
    apply WCA potential (for now)
*****************************************/

void InteractionEngine::Interact() {
  // Check if we need to update cell list
  CheckUpdate();
  // Loop through and calculate interactions
#ifdef ENABLE_OPENMP
  CalculateInteractionsMP();
#else
  CalculateInteractions();
#endif
  // Apply forces, torques, and potentials in serial
  ApplyInteractions();
  if (overlap_)
    error_exit("ERROR: Overlap of elements detected.\n");
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

int InteractionEngine::CountSpecies() {
  int obj_count = 0;
  for (auto spec=species_->begin(); spec!=species_->end(); ++spec)
    obj_count += (*spec)->GetCount();
  return obj_count;
}

void InteractionEngine::CheckUpdate() {
  // First check to see if any objects were added to the system
  int obj_count = CountSpecies();
  if (obj_count != n_objs_) {
    // reset periodic update count and update number of objects we're tracking
    i_update_ = 0; 
    n_objs_ = obj_count;
    UpdateSimples();
    UpdateInteractions();
  }
  // Otherwise check periodic update count
  else if ((++i_update_) % n_update_ == 0)
    UpdateInteractions();
}


void InteractionEngine::ProcessInteraction(std::vector<pair_interaction>::iterator pix) {
  // Avoid certain types of interactions
  auto obj1 = pix->first.first;
  auto obj2 = pix->first.second;
  Interaction *ix = &(pix->second);
  if (obj1->GetRID() == obj2->GetRID())  return;
  if (obj1->GetCID() == obj2->GetCID())  return;
  if (obj1->GetOID() == obj2->GetOID()) 
    error_exit("ERROR! Object %d attempted self-interaction!\n", obj1->GetOID());

  //interactionmindist imd;
  MinimumDistance(obj1,obj2,ix,space_);

  // XXX Don't interact if we have an overlap. This should eventually go to a max force routine
  if (ix->dr_mag2 < 0.125*SQR(obj1->GetDiameter() + obj2->GetDiameter())) {
    overlap_ = true;
    return;
  }
  // XXX Only calculate WCA potential for now
  if (ix->dr_mag2 > potentials_.wca_.GetRCut2())  return;
  potentials_.wca_.CalcPotential(ix);
  //if (ix->dr_mag2 > potentials_.r2pot_.GetRCut2()) return;
  //potentials_.r2pot_.CalcPotential(ix);
  //if (ix->dr_mag2 > potentials_.sspot_.GetRCut2()) return;
  //potentials_.sspot_.CalcPotential(ix);
}

void InteractionEngine::CalculateInteractionsMP() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<pair_interaction>::iterator, std::vector<pair_interaction>::iterator> > chunks;
  chunks.reserve(max_threads); 
  size_t chunk_size= pair_interactions_.size() / max_threads;
  auto cur_iter = pair_interactions_.begin();
  for(int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, pair_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for 
    for(int i = 0; i < max_threads; ++i) {
      for(auto pix = chunks[i].first; pix != chunks[i].second; ++pix) {
        ProcessInteraction(pix);
        // Do torque crossproducts
        cross_product(pix->second.contact1,pix->second.force,pix->second.t1,3);
        cross_product(pix->second.contact2,pix->second.force,pix->second.t2,3);
      }
    }
  }
#endif
}

void InteractionEngine::CalculateInteractions() {
  for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
    ProcessInteraction(pix);
    if (!overlap_) {
    // Do torque crossproducts
      cross_product(pix->second.contact1,pix->second.force,pix->second.t1,3);
      cross_product(pix->second.contact2,pix->second.force,pix->second.t2,3);
    }
  }
}

void InteractionEngine::ApplyInteractions() {
  for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
    auto obj1 = pix->first.first;
    auto obj2 = pix->first.second;
    Interaction *ix = &(pix->second);
    obj1->AddForce(ix->force);
    obj2->SubForce(ix->force);
    obj1->AddTorque(ix->t1);
    obj2->SubTorque(ix->t2);
    obj1->AddPotential(ix->pote);
    obj2->AddPotential(ix->pote);
    for (int i=0; i<n_dim_; ++i)
      for (int j=0; j<n_dim_; ++j)
        stress_[n_dim_*i+j] += ix->stress[n_dim_*i+j];
  }
}

// Compute pressure tensor after virial_time_avg_ steps
void InteractionEngine::CalculatePressure() {
  double inv_V = 1.0/space_->volume;
  std::fill(space_->pressure_tensor, space_->pressure_tensor+9, 0);
  // Calculate pressure tensor from stress tensor (only physical for periodic subspace)
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      // Add particle density along principle axes
      if (i==j)
        space_->pressure_tensor[n_dim_*i+j] += n_objs_*inv_V;
      // Add time-averaged virial component
      space_->pressure_tensor[n_dim_*i+j] += stress_[n_dim_*i+j]*inv_V/virial_time_avg_;
    }
  }
  // Calculate isometric pressure
  space_->pressure = 0;
  for (int i=0; i<n_dim_; ++i)
    space_->pressure += space_->pressure_tensor[n_dim_*i+i];
  space_->pressure /= n_dim_;
  // Reset local stress tensor
  std::fill(stress_,stress_+9,0);
}

bool InteractionEngine::CheckOverlap() {
  overlap_ = false;
  UpdateSimples();
  UpdateInteractions();
  CalculateInteractions();
  return overlap_;
}
