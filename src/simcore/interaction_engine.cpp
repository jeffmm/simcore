
#include "interaction_engine.h"

void InteractionEngine::Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species, 
              space_struct *space, int *i_step, bool processing) {
  // Set up pointer structures
  params_ = params;
  species_ = species;
  space_ = space;
  processing_ = processing;

  // Initialize owned structures
  no_init_ = false;
  i_step_ = i_step;
  static_pnumber_ = params_->static_particle_number;
  n_dim_ = params_->n_dim;
  n_periodic_ = params_->n_periodic;
  n_update_ = params_->n_update_cells;
  n_thermo_ = params_->n_thermo;
  no_interactions_ = !(params_->interaction_flag);
  i_update_ = -1;
  n_objs_ = -1;
  //clist_.Init(n_dim_, n_periodic_, params_->cell_length, space_->radius);
  bool local_order = (params_->local_order_analysis || params_->polar_order_analysis || params_->overlap_analysis || params_->density_analysis);
  if (processing && local_order) {
    params_->cell_length = 0.5*params_->local_order_width;
  }
  ptracker_.Init(params, &interactors_, &nlist_);
  // Update dr distance should be half the cell length, and we are comparing the squares of the trajectory distances
  dr_update_ = 0.25*ptracker_.GetCellLength()*ptracker_.GetCellLength();
  mindist_.Init(space, 2.0*dr_update_);
  potentials_.InitPotentials(params_);
  if (local_order && processing) {
    struct_analysis_.Init(params, i_step);
  }
  xlink_.Init(params_, space_, &mindist_, &ix_objects_);
}

/****************************************
  INTERACT: Loop through interactions and
    apply WCA potential (for now)
*****************************************/

void InteractionEngine::Interact() {
  n_interactions_ = 0;
  // First check if we need to interact
  if (no_interactions_ && space_->type == +boundary_type::none) return;
  // Check if we need to update cell list
  xlink_.UpdateCrosslinks();
  CheckUpdate();
  /* Update anchors and crosslinks */
  // Loop through and calculate interactions
  if (! no_interactions_ ) {
    CalculatePairInteractions();
  }
  CalculateBoundaryInteractions();
  // Apply forces, torques, and potentials in serial
  if (! no_interactions_ ) {
    ApplyPairInteractions();
  }
  ApplyBoundaryInteractions();


  if (params_->in_out_flag) {
    if (n_interactions_ == 0 && in_out_flag_) {
      early_exit = true;
    }
    else if (n_interactions_ > 0 && !in_out_flag_) {
      in_out_flag_ = true;
    }
  }
}

void InteractionEngine::ForceUpdate() {
  UpdateInteractors();
  UpdateInteractions();
}

void InteractionEngine::UpdateInteractors() {
  ix_objects_.clear();
  interactors_.clear();
  std::vector<Object*> ixs;
  for (auto spec_it= species_->begin(); spec_it!=species_->end(); ++spec_it) {
    (*spec_it)->GetInteractors(&ixs);
    ix_objects_.insert(ix_objects_.end(), ixs.begin(), ixs.end());
  }
  interactors_.insert(interactors_.end(), ix_objects_.begin(), ix_objects_.end());
  // Add crosslinks as interactors
  std::vector<Object*> xlinks;
  xlink_.GetInteractors(&xlinks);
  interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
}

void InteractionEngine::UpdateInteractions() {
  if (no_interactions_) return;
  ptracker_.AssignCells();
  ptracker_.CreatePairsCellList();
  pair_interactions_.clear();
  for (auto pair=nlist_.begin(); pair!=nlist_.end(); ++pair) {
    Interaction ix;
    pair_interaction pix(std::make_pair(&(*(interactors_[pair->first])), &(*(interactors_[pair->second]))), ix);
    pair_interactions_.push_back(pix);
  }
  boundary_interactions_.clear();
  for (auto ixor=interactors_.begin(); ixor!=interactors_.end(); ++ixor) {
    Interaction ix;
    if (mindist_.CheckBoundaryInteraction(*ixor, &ix)) {
      boundary_interactions_.push_back(std::make_pair(*ixor, ix));
    }
  }

}

int InteractionEngine::CountSpecies() {
  int obj_count = 0;
  for (auto spec=species_->begin(); spec!=species_->end(); ++spec) {
    obj_count += (*spec)->GetCount();
  }
  return obj_count;
}

bool InteractionEngine::CountAndUpdate() {
  if (processing_) return;
  int obj_count = CountSpecies();
  if (obj_count != n_objs_) {
    // reset update count and update number of objects to track
    i_update_ = 0; 
    n_objs_ = obj_count;
    ForceUpdate();
    xlink_.UpdateObjsVolume();
    return true;
  }
  return false;
}

void InteractionEngine::CheckUpdate() {
  /* First check to see if any objects were added to the system;
     If static_pnumber_ is flagged, we know that particle numbers
     never change, so we don't bother counting particles and move on */
  if (!static_pnumber_) {
    if (CountAndUpdate()) return;
  }

  /* If n_update_ <=0, we update nearest neighbors if any particle 
     has moved a distance further than dr_update_ */
  if (n_update_ <=0 && GetDrMax() > dr_update_) {
    i_update_ = 0; 
    UpdateInteractions();
    ZeroDrTot();
  }

  /* If n_update_ > 0, we use that number as an update frequency,
     and update nearest neighbors every n_update_ steps */
  else if (n_update_ > 0 && (++i_update_) % n_update_ == 0) {
    UpdateInteractions();
  }

  /* If we need to update crosslinks */
  if (xlink_.CheckUpdate()) {
    /*TODO This might be too frequent: can I do this better?*/
    ForceUpdate();
    /*TODO Figure out why the below code leads to memory leaks */
    //interactors_.clear();
    //interactors_.insert(interactors_.end(), ix_objects_.begin(), ix_objects_.end());
    /* Add crosslinks as interactors */
    //std::vector<Object*> xlinks;
    //xlink_.GetInteractors(&xlinks);
    //interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
  }
}

/* Returns maximum distance traveled by any particle of any species in the system */
double InteractionEngine::GetDrMax() {
  double max_dr = 0;
  for (auto spec=species_->begin(); spec!=species_->end(); ++spec) {
    double dr = (*spec)->GetDrMax();
    if (dr>max_dr) {
      max_dr = dr;
    }
  }
  return max_dr;
}

/* Resets the origin of each particle's tracked trajectory, for tracking
   total distance travelled by particle*/
void InteractionEngine::ZeroDrTot() {
  for (auto spec=species_->begin(); spec!=species_->end(); ++spec) {
    (*spec)->ZeroDrTot();
  }
}

void InteractionEngine::ProcessPairInteraction(std::vector<pair_interaction>::iterator pix) {
  // Avoid certain types of interactions
  auto obj1 = pix->first.first;
  auto obj2 = pix->first.second;
  Interaction *ix = &(pix->second);
  // Rigid objects don't self interact
  //if (obj1->GetRID() == obj2->GetRID())  return;
  // Composite objects do self interact if they want to
  // Check to make sure we aren't self-interacting at the simple object level
  if (obj1->GetOID() == obj2->GetOID()) {
    //error_exit("Object %d attempted self-interaction!", obj1->GetOID());
  }
  // If we are disallowing like-like interactions, then similar species do not interact...
  // ...so do not interact if same species
  if (!params_->like_like_interactions && obj1->GetSID() == obj2->GetSID()) {
    return;
  }
  // Check that objects are not both motors, which do not interact
  if (obj1->GetSID() == +species_id::crosslink && obj2->GetSID() == +species_id::crosslink) {
    return;
  }

  // Check that object 1 is part of a mesh, in which case...
  // ...check that object 1 is of the same mesh of object 2, in which case...
  if (obj1->GetMeshID() > 0 && obj1->GetMeshID() == obj2->GetMeshID()) {
    // ..check if object 1 or object 2 are crosslinks, in which case: do not interact
    if (obj1->GetSID() == +species_id::crosslink || obj2->GetSID() == +species_id::crosslink) {
      return;
    }

    // ...check if object 2 is a neighbor of object 1, in which case: do not interact
    if (obj1->HasNeighbor(obj2->GetOID())) {
      return;
    }
  }

  if (params_->filament.spiral_flag == 1 && obj1->GetMeshID() != obj2->GetMeshID()) {
    return;
  }

  /* If one object is a crosslink, add object to crosslink neighbor list */
  if (obj1->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToXlink(obj1, obj2);
    return;
  }
  else if (obj2->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToXlink(obj2, obj1);
    return;
  }

  // We have an interaction:
  if ( obj1->GetMeshID() != obj2->GetMeshID() ) {
    n_interactions_++;
  }
  if (obj1->GetSID() == +species_id::crosslink || obj2->GetSID() == +species_id::crosslink) {
    //printf("Crosslink interacting\n");
    //printf("   with: %s and %s \n", obj1->GetSID()._to_string(), obj2->GetSID()._to_string());
  }

  mindist_.ObjectObject(obj1, obj2, ix);

  // XXX Don't interact if we have an overlap. This should eventually go to a max force routine
  if (ix->dr_mag2 < 0.25*SQR(obj1->GetDiameter() + obj2->GetDiameter())) {
    overlap_ = true;
  }
  /* Check to see if particles are not close enough to interact */
  if (ix->dr_mag2 > potentials_.GetRCut2())  return;
  /* Calculates forces from the potential defined during initialization */
  potentials_.CalcPotential(ix);
}

void InteractionEngine::ProcessBoundaryInteraction(std::vector<boundary_interaction>::iterator bix) {
  mindist_.BoundaryCondition(bix);
  Interaction * ix = &(bix->second);
  // XXX Don't interact if we have an overlap. This should eventually go to a max force routine
  if (ix->dr_mag2 < 0.25*SQR(bix->first->GetDiameter())) {
    overlap_ = true;
  }
  if (ix->dr_mag2 > potentials_.GetRCut2()) return;
  potentials_.CalcPotential(ix);
}

void InteractionEngine::CalculateBoundaryInteractions() {
  if (space_->type == +boundary_type::none) {
    return;
  }
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<boundary_interaction>::iterator, std::vector<boundary_interaction>::iterator> > chunks;
  chunks.reserve(max_threads); 
  size_t chunk_size= boundary_interactions_.size() / max_threads;
  auto cur_iter = boundary_interactions_.begin();
  for(int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, boundary_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for 
    for(int i = 0; i < max_threads; ++i) {
      for(auto bix = chunks[i].first; bix != chunks[i].second; ++bix) {
        ProcessBoundaryInteraction(bix);
        // Do torque crossproducts
        cross_product(bix->second.contact1, bix->second.force, bix->second.t1, 3);
      }
    }
  }
#else 
  for(auto bix = boundary_interactions_.begin(); bix != boundary_interactions_.end(); ++bix) {
    ProcessBoundaryInteraction(bix);
    // Do torque crossproducts
    cross_product(bix->second.contact1, bix->second.force, bix->second.t1, 3);
  }
#endif
}

void InteractionEngine::CalculatePairInteractions() {
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
        ProcessPairInteraction(pix);
        // Do torque crossproducts
        cross_product(pix->second.contact1, pix->second.force, pix->second.t1, 3);
        cross_product(pix->second.contact2, pix->second.force, pix->second.t2, 3);
      }
    }
  }
#else
  for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
    ProcessPairInteraction(pix);
    // Do torque crossproducts
    cross_product(pix->second.contact1, pix->second.force, pix->second.t1, 3);
    cross_product(pix->second.contact2, pix->second.force, pix->second.t2, 3);
  }
#endif
}

void InteractionEngine::ApplyPairInteractions() {
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
    for (int i=0; i<n_dim_; ++i) {
      for (int j=0; j<n_dim_; ++j) {
        stress_[n_dim_*i+j] += ix->stress[n_dim_*i+j];
      }
    }
  }
}

void InteractionEngine::ApplyBoundaryInteractions() {
  for(auto bix = boundary_interactions_.begin(); bix != boundary_interactions_.end(); ++bix) {
    auto obj1 = bix->first;
    Interaction *ix = &(bix->second);
    obj1->AddForce(ix->force);
    obj1->AddTorque(ix->t1);
    obj1->AddPotential(ix->pote);
    for (int i=0; i<n_dim_; ++i) {
      for (int j=0; j<n_dim_; ++j) {
        stress_[n_dim_*i+j] += ix->stress[n_dim_*i+j];
      }
    }
  }
}

// Compute pressure tensor after n_thermo_ steps
void InteractionEngine::CalculatePressure() {
  double inv_V = 1.0/space_->volume;
  std::fill(space_->pressure_tensor, space_->pressure_tensor+9, 0);
  // Calculate pressure tensor from stress tensor (only physical for periodic subspace)
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      // Add particle density along principle axes
      if (i==j) {
        space_->pressure_tensor[n_dim_*i+j] += n_objs_*inv_V;
      }
      // Add time-averaged virial component
      space_->pressure_tensor[n_dim_*i+j] += stress_[n_dim_*i+j]*inv_V/n_thermo_;
    }
  }
  // Calculate isometric pressure
  space_->pressure = 0;
  for (int i=0; i<n_dim_; ++i) {
    space_->pressure += space_->pressure_tensor[n_dim_*i+i];
  }
  space_->pressure /= n_dim_;
  // Reset local stress tensor
  std::fill(stress_, stress_+9, 0);
}

bool InteractionEngine::CheckOverlap(std::vector<Object*> ixs) {
  overlap_ = false;
  /* Only consider objects (not crosslinks) for overlaps */
  int n_interactors = ix_objects_.size();
  ptracker_.CreatePartialPairsCellList(ixs, n_interactors);
  pair_interactions_.clear();
  for (auto pair=nlist_.begin(); pair!=nlist_.end(); ++pair) {
    Object *o1, *o2;
    Interaction ix;
    if (pair->first > n_interactors-1 && pair->second < n_interactors) {
      o1 = ixs[pair->first - n_interactors];
      o2 = ix_objects_[pair->second];
    }
    else if (pair->second > n_interactors-1 && pair->first < n_interactors) {
      o1 = ixs[pair->second - n_interactors];
      o2 = ix_objects_[pair->first];
    }
    else if (pair->second > n_interactors-1 && pair->first > n_interactors-1) {
      o1 = ixs[pair->first - n_interactors];
      o2 = ixs[pair->second - n_interactors];
    }
    else {
      o1 = ix_objects_[pair->first];
      o2 = ix_objects_[pair->second];
    }
    pair_interaction pix(std::make_pair(o1, o2), ix);
    pair_interactions_.push_back(pix);
  }
  CalculatePairInteractions();
  return overlap_;
}

bool InteractionEngine::CheckBoundaryConditions(std::vector<Object*> ixs) {
  bool outside_boundary = false;
  for (auto ixor=ixs.begin(); ixor!=ixs.end(); ++ixor) {
    if (mindist_.CheckOutsideBoundary(*ixor)) {
      outside_boundary = true;
    }
  }
  return outside_boundary;
}

/* Here, I want to calculate P(r, phi), which tells me the probability
   of finding an object at a position (r, phi) in its reference frame. */
void InteractionEngine::StructureAnalysis() {
  ForceUpdate();
  CalculateStructure();
}

void InteractionEngine::CalculateStructure() {
  /* Only consider objects (not crosslinks) for structure analysis */
  if (struct_analysis_.GetNumObjs() == 0) {
    int nobj=CountSpecies();
    struct_analysis_.SetNumObjs(nobj);
  }
  struct_analysis_.IncrementCount();
  if (params_->polar_order_analysis) {
    for (auto it=ix_objects_.begin(); it!=ix_objects_.end(); ++it) {
      (*it)->ZeroPolarOrder();
    }
  }
  if (!no_interactions_) {
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
          if (params_->polar_order_analysis || params_->overlap_analysis) {
            mindist_.ObjectObject(pix->first.first, pix->first.second, &(pix->second));
          }
          struct_analysis_.CalculateStructurePair(pix);
        }
      }
    }
#else
    for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
      if (params_->polar_order_analysis || params_->overlap_analysis) {
        mindist_.ObjectObject(pix->first.first, pix->first.second, &(pix->second));
      }
      struct_analysis_.CalculateStructurePair(pix);
    }
#endif
  }
  //if (params_->overlap_analysis) {
    //for(auto pix = pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
      //struct_analysis_.CountOverlap(pix);
    //}
  //}
  struct_analysis_.AverageStructure();
  if (params_->polar_order_analysis) {
    for (auto pix=pair_interactions_.begin(); pix != pair_interactions_.end(); ++pix) {
      auto obj1 = pix->first.first;
      auto obj2 = pix->first.second;
      obj1->AddPolarOrder(pix->second.polar_order);
      obj2->AddPolarOrder(pix->second.polar_order);
      obj1->AddContactNumber(pix->second.contact_number);
      obj2->AddContactNumber(pix->second.contact_number);
    }
    for (auto it=ix_objects_.begin(); it!=ix_objects_.end(); ++it) {
      (*it)->CalcPolarOrder();
    }
  }
  if (params_->density_analysis) {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    std::vector<std::pair<std::vector<Object*>::iterator, std::vector<Object*>::iterator> > chunks;
    chunks.reserve(max_threads); 
    size_t chunk_size= ix_objects_.size() / max_threads;
    auto cur_iter = ix_objects_.begin();
    for(int i = 0; i < max_threads - 1; ++i) {
      auto last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, ix_objects_.end()));
#pragma omp parallel shared(chunks)
    {
#pragma omp for 
      for(int i = 0; i < max_threads; ++i) {
        for(auto it = chunks[i].first; it != chunks[i].second; ++it) {
          struct_analysis_.BinDensity(*it);
        }
      }
    }
#else
    for (auto it=ix_objects_.begin(); it!=ix_objects_.end(); ++it) {
      struct_analysis_.BinDensity(*it);
    }
#endif
  }
}

void InteractionEngine::Clear() {
  if (no_init_) return;
  ptracker_.Clear();
  bool local_order = (params_->local_order_analysis || 
      params_->polar_order_analysis || 
      params_->overlap_analysis || 
      params_->density_analysis);
  if (local_order && processing_) {
    struct_analysis_.Clear();
  }
  xlink_.Clear();
}

/* Only used during species insertion */
void InteractionEngine::Reset() {
  pair_interactions_.clear();
  ix_objects_.clear();
  interactors_.clear();
  ptracker_.ClearCells();
}

/* Only used during species insertion */
void InteractionEngine::AddInteractors(std::vector<Object*> ixs) {
  ptracker_.AddToCellList(ixs, ix_objects_.size());
  ix_objects_.insert(ix_objects_.end(), ixs.begin(), ixs.end());
}

void InteractionEngine::DrawInteractions(std::vector<graph_struct*> * graph_array) {
  xlink_.Draw(graph_array);
}

void InteractionEngine::WriteOutputs() {
  xlink_.WriteOutputs();
}

void InteractionEngine::InitOutputs(bool reading_inputs, bool reduce_flag, bool with_reloads) {
  xlink_.InitOutputs(reading_inputs, reduce_flag, with_reloads);
  if (params_->load_checkpoint) {
    /* TODO Run minimum distance to determine which anchor belongs to which bond */ 
  }
}

void InteractionEngine::ReadInputs() {
  xlink_.ReadInputs();
}
