
#include "interaction_engine.hpp"

void InteractionEngine::Init(system_parameters *params,
                             std::vector<SpeciesBase *> *species,
                             space_struct *space, int *i_step,
                             bool processing) {
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
  int n_cells_1d =
      (int)floor(2 * params_->system_radius / params_->cell_length);
#ifdef TRACE
  if (n_cells_1d > 20) {
    Logger::Warning("Simulation run in trace mode with a large number of "
                    "cells in cell list (%d).",
                    n_cells_1d);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif
  double cell_length = (double)2 * params_->system_radius / n_cells_1d;
  clist_.Init(n_cells_1d, cell_length, params_->n_dim, params_->n_periodic);
  bool local_order =
      (params_->local_order_analysis || params_->polar_order_analysis ||
       params_->overlap_analysis || params_->density_analysis);
  if (processing && local_order) {
    params_->cell_length = 0.5 * params_->local_order_width;
  }
  // Update dr distance should be half the cell length, and we are comparing the
  // squares of the trajectory distances
  dr_update_ = 0.25 * cell_length * cell_length;
  mindist_.Init(space, 2.0 * dr_update_);
  potentials_.InitPotentials(params_);
  if (local_order && processing) {
    struct_analysis_.Init(params, i_step);
  }
  xlink_.Init(params_, space_, &mindist_, &ix_objects_);
  no_boundaries_ = false;
  if (space_->type == +boundary_type::none)
    no_boundaries_ = true;
}

/****************************************
  INTERACT: Loop through interactions and
    apply WCA potential (for now)
*****************************************/

void InteractionEngine::Interact() {
  n_interactions_ = 0;
  // First check if we need to interact
  if (no_interactions_ && no_boundaries_)
    return;
  // Check if we need to update objects in cell list
  CheckUpdateObjects();
  // Update crosslinks
  xlink_.UpdateCrosslinks();
  // Check if we need to update crosslink interactors
  CheckUpdateXlinks();
  CheckUpdateInteractions();

  /* Update anchors and crosslinks */
  // Loop through and calculate interactions
  if (!no_interactions_) {
    CalculatePairInteractions();
  }
  CalculateBoundaryInteractions();
  // Apply forces, torques, and potentials in serial
  if (!no_interactions_) {
    ApplyPairInteractions();
  }
  ApplyBoundaryInteractions();

  if (params_->in_out_flag) {
    if (n_interactions_ == 0 && in_out_flag_) {
      early_exit = true;
    } else if (n_interactions_ > 0 && !in_out_flag_) {
      in_out_flag_ = true;
    }
  }
}

void InteractionEngine::CheckUpdateXlinks() {
  /* If we need to update crosslinks */
  if (xlink_.CheckUpdate()) {
    interactors_.clear();
    interactors_.insert(interactors_.end(), ix_objects_.begin(),
                        ix_objects_.end());
    // Add crosslinks as interactors
    std::vector<Object *> xlinks;
    xlink_.GetInteractors(xlinks);
    interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
    UpdateInteractions();
    i_update_ = 0;
    if (n_update_ <= 0) {
      ZeroDrTot();
    }
    return;
  }
}
void InteractionEngine::ForceUpdate() {
  UpdateInteractors();
  UpdateInteractions();
}

void InteractionEngine::UpdateInteractors() {
  ix_objects_.clear();
  interactors_.clear();
  for (auto spec_it = species_->begin(); spec_it != species_->end();
       ++spec_it) {
    (*spec_it)->GetInteractors(ix_objects_);
  }
  // Add crosslinks as interactors
  interactors_.insert(interactors_.end(), ix_objects_.begin(),
                      ix_objects_.end());
  std::vector<Object *> xlinks;
  xlink_.GetInteractors(xlinks);
  interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
}

/* Checks whether or not the given anchor is supposed to be attached to the
   given bond by checking mesh_id of both the bond and the anchor. If they
   match, the anchor is attached to the mesh using the anchor mesh_lambda */
bool InteractionEngine::CheckBondAnchorPair(Object *anchor, Object *bond) {
  // Check that the bond and anchor share a mesh_id
  if (anchor->GetMeshID() == bond->GetMeshID()) {
    Anchor *a = dynamic_cast<Anchor *>(anchor);
    if (a == nullptr) {
      Logger::Error("Object pointer was unsuccessfully dynamically cast to an "
                    "Anchor pointer in CheckBondAnchorPair!");
    }
    if (bond->GetType() != +obj_type::bond) {
      Logger::Error(
          "CheckBondAnchorPair expected Bond object pointer, but object"
          " pointer does not have bond obj_type!");
    }
    // Check that the anchor isn't already attached
    if (a->GetBondLambda() < 0) {
      a->AttachObjMeshLambda(bond, a->GetMeshLambda());
      return true;
    }
  }
  return false;
}

void InteractionEngine::PairBondCrosslinks() {
  Logger::Trace("Pairing bound crosslinks and objects");
  ix_objects_.clear();
  interactors_.clear();
  // First get object interactors (non-crosslinks)
  for (auto spec_it = species_->begin(); spec_it != species_->end();
       ++spec_it) {
    (*spec_it)->GetInteractors(ix_objects_);
  }
  interactors_.insert(interactors_.end(), ix_objects_.begin(),
                      ix_objects_.end());

  // Add anchors as interactors
  std::vector<Object *> anchors;
  xlink_.GetAnchorInteractors(anchors);
  interactors_.insert(interactors_.end(), anchors.begin(), anchors.end());
  pair_interactions_.clear();
  clist_.RenewObjectsCells(interactors_);
  clist_.MakePairs(pair_interactions_);
  int n_anchors_attached = 0;
  // for (auto ix = anchors.begin(); ix != anchors.end(); ++ix) {
  // for (auto jx = ix_objects_.begin(); jx != ix_objects_.end(); ++jx) {
  for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
       ++ix) {
    Object *obj1 = ix->obj1;
    Object *obj2 = ix->obj2;
    if (obj1->GetSID() == +species_id::crosslink &&
        obj2->GetType() == +obj_type::bond) {
      if (CheckBondAnchorPair(obj1, obj2))
        // if (CheckBondAnchorPair(*ix, *jx))
        n_anchors_attached++;
    } else if (obj2->GetSID() == +species_id::crosslink &&
               obj1->GetType() == +obj_type::bond) {
      if (CheckBondAnchorPair(obj2, obj1)) {
        n_anchors_attached++;
      }
    }
  }
  /* Check that all anchors found their bond attachments */
  if (n_anchors_attached != anchors.size()) {
    Logger::Error("Not all anchors have found their bond! %d/%lu",
                  n_anchors_attached, anchors.size());
  }
}

void InteractionEngine::UpdateInteractions() {
  UpdatePairInteractions();
  UpdateBoundaryInteractions();
}

void InteractionEngine::UpdatePairInteractions() {
  if (no_interactions_)
    return;
  pair_interactions_.clear();
  clist_.RenewObjectsCells(interactors_);
  clist_.MakePairs(pair_interactions_);
}

void InteractionEngine::UpdateBoundaryInteractions() {
  if (no_boundaries_)
    return;
  boundary_interactions_.clear();
  for (auto ixor = interactors_.begin(); ixor != interactors_.end(); ++ixor) {
    Interaction ix(*ixor);
    if (mindist_.CheckBoundaryInteraction(ix)) {
      boundary_interactions_.push_back(ix);
    }
  }
}

int InteractionEngine::CountSpecies() {
  int obj_count = 0;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    obj_count += (*spec)->GetCount();
  }
  return obj_count;
}

const bool InteractionEngine::CheckSpeciesInteractorUpdate() const {
  bool result = false;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->CheckInteractorUpdate()) {
      result = true;
    }
  }
  return result;
}

void InteractionEngine::CheckUpdateObjects() {
  /* First check to see if any objects were added to the system;
     If static_pnumber_ is flagged, we know that particle numbers
     never change, so we don't bother counting particles and move on */
  if (processing_)
    return;
  if (static_pnumber_)
    return;
  bool ix_update = CheckSpeciesInteractorUpdate();
  int obj_count = CountSpecies();
  if (obj_count != n_objs_ || ix_update) {
    // reset update count and update number of objects to track
    i_update_ = 0;
    n_objs_ = obj_count;
    ForceUpdate();
    xlink_.UpdateObjsVolume();
    // printf("Forcing interactor update, i_step: %d\n", params_->i_step);
  }
}

void InteractionEngine::ResetCellList() { clist_.ResetNeighbors(); }

void InteractionEngine::CheckUpdateInteractions() {
  /* If n_update_ <=0, we update nearest neighbors if any particle
     has moved a distance further than dr_update_ */
  if (n_update_ <= 0 && GetDrMax() > dr_update_) {
    i_update_ = 0;
    UpdateInteractions();
    ZeroDrTot();
  }
  /* If n_update_ > 0, we use that number as an update frequency,
     and update nearest neighbors every n_update_ steps */
  else if (n_update_ > 0 && (++i_update_) % n_update_ == 0) {
    UpdateInteractions();
  }
}

/* Returns maximum distance traveled by any particle of any species in the
 * system */
double InteractionEngine::GetDrMax() {
  double max_dr = 0;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    double dr = (*spec)->GetDrMax();
    if (dr > max_dr) {
      max_dr = dr;
    }
  }
  return max_dr;
}

/* Resets the origin of each particle's tracked trajectory, for tracking
   total distance travelled by particle*/
void InteractionEngine::ZeroDrTot() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->ZeroDrTot();
  }
}

void InteractionEngine::ProcessPairInteraction(ix_iterator ix) {
  // Avoid certain types of interactions
  Object *obj1 = ix->obj1;
  Object *obj2 = ix->obj2;
  Logger::Trace("Processing interaction between %d and %d", obj1->GetOID(),
                obj2->GetOID());
  // Rigid objects don't self interact
  // if (obj1->GetRID() == obj2->GetRID())  return;
  // Composite objects do self interact if they want to
  // Check to make sure we aren't self-interacting at the simple object level
  if (obj1->GetOID() == obj2->GetOID()) {
    // Logger::Error("Object %d attempted self-interaction!", obj1->GetOID());
  }
  // If we are disallowing like-like interactions, then similar species do not
  // interact...
  // ...so do not interact if same species
  if (!params_->like_like_interactions && obj1->GetSID() == obj2->GetSID()) {
    return;
  }
  // Check that objects are not both motors, which do not interact
  if (obj1->GetSID() == +species_id::crosslink &&
      obj2->GetSID() == +species_id::crosslink) {
    return;
  }

  // Check that object 1 is part of a mesh, in which case...
  // ...check that object 1 is of the same mesh of object 2, in which case...
  if (obj1->GetMeshID() > 0 && obj1->GetMeshID() == obj2->GetMeshID()) {
    // ..check if object 1 or object 2 are crosslinks, in which case: do not
    // interact
    if (obj1->GetSID() == +species_id::crosslink ||
        obj2->GetSID() == +species_id::crosslink) {
      return;
    }

    // ...check if object 2 is a neighbor of object 1, in which case: do not
    // interact
    if (obj1->HasNeighbor(obj2->GetOID())) {
      return;
    }
  }

  /* If one object is a crosslink, add object to crosslink neighbor list */
  if (obj1->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToAnchor(obj1, obj2);
    return;
  } else if (obj2->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToAnchor(obj2, obj1);
    return;
  }

  // We have an interaction:
  if (obj1->GetMeshID() != obj2->GetMeshID()) {
    n_interactions_++;
  }
  if (obj1->GetSID() == +species_id::crosslink ||
      obj2->GetSID() == +species_id::crosslink) {
    // printf("Crosslink interacting\n");
    // printf("   with: %s and %s \n", obj1->GetSID()._to_string(),
    // obj2->GetSID()._to_string());
  }

  mindist_.ObjectObject(*ix);

  // XXX Don't interact if we have an overlap. This should eventually go to a
  // max force routine
  if (ix->dr_mag2 < 0.25 * SQR(obj1->GetDiameter() + obj2->GetDiameter())) {
    overlap_ = true;
  }
  /* Check to see if particles are not close enough to interact */
  if (ix->dr_mag2 > potentials_.GetRCut2())
    return;
  /* Calculates forces from the potential defined during initialization */
  potentials_.CalcPotential(*ix);
}

void InteractionEngine::ProcessBoundaryInteraction(ix_iterator ix) {
  mindist_.CheckBoundaryInteraction(*ix);
  // XXX Don't interact if we have an overlap. This should eventually go to a
  // max force routine
  if (ix->dr_mag2 < 0.25 * SQR(ix->obj1->GetDiameter())) {
    overlap_ = true;
  }
  if (ix->dr_mag2 > potentials_.GetRCut2())
    return;
  potentials_.CalcPotential(*ix);
}

void InteractionEngine::CalculateBoundaryInteractions() {
  if (space_->type == +boundary_type::none) {
    return;
  }
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<ix_iterator, ix_iterator>> chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = boundary_interactions_.size() / max_threads;
  auto cur_iter = boundary_interactions_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, boundary_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto ix = chunks[i].first; ix != chunks[i].second; ++ix) {
        ProcessBoundaryInteraction(ix);
        // Do torque crossproducts
        cross_product(ix->contact1, ix->force, ix->t1, 3);
      }
    }
  }
#else
  for (auto ix = boundary_interactions_.begin();
       ix != boundary_interactions_.end(); ++ix) {
    ProcessBoundaryInteraction(ix);
    // Do torque crossproducts
    cross_product(ix->contact1, ix->force, ix->t1, 3);
  }
#endif
}

void InteractionEngine::CalculatePairInteractions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<ix_iterator, ix_iterator>> chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = pair_interactions_.size() / max_threads;
  auto cur_iter = pair_interactions_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, pair_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto ix = chunks[i].first; ix != chunks[i].second; ++ix) {
        ProcessPairInteraction(ix);
        // Do torque crossproducts
        cross_product(ix->contact1, ix->force, ix->t1, 3);
        cross_product(ix->contact2, ix->force, ix->t2, 3);
      }
    }
  }
#else
  for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
       ++ix) {
    ProcessPairInteraction(ix);
    // Do torque crossproducts
    cross_product(ix->contact1, ix->force, ix->t1, 3);
    cross_product(ix->contact2, ix->force, ix->t2, 3);
  }
#endif
}

void InteractionEngine::ApplyPairInteractions() {
  for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
       ++ix) {
    Object *obj1 = ix->obj1;
    Object *obj2 = ix->obj2;
    obj1->AddForce(ix->force);
    obj2->SubForce(ix->force);
    obj1->AddTorque(ix->t1);
    obj2->SubTorque(ix->t2);
    obj1->AddPotential(ix->pote);
    obj2->AddPotential(ix->pote);
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        stress_[n_dim_ * i + j] += ix->stress[n_dim_ * i + j];
      }
    }
  }
}

void InteractionEngine::ApplyBoundaryInteractions() {
  for (auto ix = boundary_interactions_.begin();
       ix != boundary_interactions_.end(); ++ix) {
    Object *obj1 = ix->obj1;
    obj1->AddForce(ix->force);
    obj1->AddTorque(ix->t1);
    obj1->AddPotential(ix->pote);
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        stress_[n_dim_ * i + j] += ix->stress[n_dim_ * i + j];
      }
    }
  }
}

// Compute pressure tensor after n_thermo_ steps
void InteractionEngine::CalculatePressure() {
  double inv_V = 1.0 / space_->volume;
  std::fill(space_->pressure_tensor, space_->pressure_tensor + 9, 0);
  // Calculate pressure tensor from stress tensor (only physical for periodic
  // subspace)
  for (int i = 0; i < n_dim_; ++i) {
    for (int j = 0; j < n_dim_; ++j) {
      // Add particle density along principle axes
      if (i == j) {
        space_->pressure_tensor[n_dim_ * i + j] += n_objs_ * inv_V;
      }
      // Add time-averaged virial component
      space_->pressure_tensor[n_dim_ * i + j] +=
          stress_[n_dim_ * i + j] * inv_V / n_thermo_;
    }
  }
  // Calculate isometric pressure
  space_->pressure = 0;
  for (int i = 0; i < n_dim_; ++i) {
    space_->pressure += space_->pressure_tensor[n_dim_ * i + i];
  }
  space_->pressure /= n_dim_;
  // Reset local stress tensor
  std::fill(stress_, stress_ + 9, 0);
}

bool InteractionEngine::CheckOverlap(std::vector<Object *> &ixors) {
  overlap_ = false;
  /* Only consider objects (not crosslinks) for overlaps */
  for (auto ixor = ixors.begin(); ixor != ixors.end(); ++ixor) {
    pair_interactions_.clear();
    clist_.PairSingleObject(**ixor, pair_interactions_);
    CalculatePairInteractions();
    if (overlap_)
      break;
  }
  return overlap_;
}

bool InteractionEngine::CheckBoundaryConditions(std::vector<Object *> &ixors) {
  bool outside_boundary = false;
  for (auto ixor = ixors.begin(); ixor != ixors.end(); ++ixor) {
    if (mindist_.CheckOutsideBoundary(**ixor)) {
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
    int nobj = CountSpecies();
    struct_analysis_.SetNumObjs(nobj);
  }
  struct_analysis_.IncrementCount();
  if (params_->polar_order_analysis) {
    for (auto it = ix_objects_.begin(); it != ix_objects_.end(); ++it) {
      (*it)->ZeroPolarOrder();
    }
  }
  if (!no_interactions_) {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    std::vector<std::pair<ix_iterator, ix_iterator>> chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = pair_interactions_.size() / max_threads;
    auto cur_iter = pair_interactions_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      auto last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, pair_interactions_.end()));
#pragma omp parallel shared(chunks)
    {
#pragma omp for
      for (int i = 0; i < max_threads; ++i) {
        for (auto ix = chunks[i].first; ix != chunks[i].second; ++ix) {
          if (params_->polar_order_analysis || params_->overlap_analysis) {
            mindist_.ObjectObject(*ix);
          }
          struct_analysis_.CalculateStructurePair(ix);
        }
      }
    }
#else
    for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
         ++ix) {
      if (params_->polar_order_analysis || params_->overlap_analysis) {
        mindist_.ObjectObject(*ix);
      }
      struct_analysis_.CalculateStructurePair(ix);
    }
#endif
  }
  // if (params_->overlap_analysis) {
  // for(auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
  // ++ix) { struct_analysis_.CountOverlap(ix);
  //}
  //}
  struct_analysis_.AverageStructure();
  if (params_->polar_order_analysis) {
    for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
         ++ix) {
      Object *obj1 = ix->obj1;
      Object *obj2 = ix->obj2;
      obj1->AddPolarOrder(ix->polar_order);
      obj2->AddPolarOrder(ix->polar_order);
      obj1->AddContactNumber(ix->contact_number);
      obj2->AddContactNumber(ix->contact_number);
    }
    for (auto it = ix_objects_.begin(); it != ix_objects_.end(); ++it) {
      (*it)->CalcPolarOrder();
    }
  }
  if (params_->density_analysis) {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    std::vector<std::pair<std::vector<Object *>::iterator,
                          std::vector<Object *>::iterator>>
        chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = ix_objects_.size() / max_threads;
    auto cur_iter = ix_objects_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      auto last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, ix_objects_.end()));
#pragma omp parallel shared(chunks)
    {
#pragma omp for
      for (int i = 0; i < max_threads; ++i) {
        for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
          struct_analysis_.BinDensity(*it);
        }
      }
    }
#else
    for (auto it = ix_objects_.begin(); it != ix_objects_.end(); ++it) {
      struct_analysis_.BinDensity(*it);
    }
#endif
  }
}

void InteractionEngine::Clear() {
  if (no_init_)
    return;
  clist_.Clear();
  bool local_order =
      (params_->local_order_analysis || params_->polar_order_analysis ||
       params_->overlap_analysis || params_->density_analysis);
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
  clist_.ClearCellObjects();
}

/* Only used during species insertion */
void InteractionEngine::AddInteractors(std::vector<Object *> &ixs) {
  clist_.AssignObjectsCells(ixs);
  ix_objects_.insert(ix_objects_.end(), ixs.begin(), ixs.end());
}

void InteractionEngine::DrawInteractions(
    std::vector<graph_struct *> &graph_array) {
  xlink_.Draw(graph_array);
}

void InteractionEngine::WriteOutputs() { xlink_.WriteOutputs(); }

void InteractionEngine::InitOutputs(bool reading_inputs, bool reduce_flag,
                                    bool with_reloads) {
  xlink_.InitOutputs(reading_inputs, reduce_flag, with_reloads);
  if (params_->load_checkpoint) {
    PairBondCrosslinks();
    ForceUpdate();
  }
}

void InteractionEngine::ReadInputs() { xlink_.ReadInputs(); }

void InteractionEngine::InitCrosslinkSpecies(sid_label &slab,
                                             ParamsParser &parser) {
  xlink_.InitSpecies(slab, parser);
}
