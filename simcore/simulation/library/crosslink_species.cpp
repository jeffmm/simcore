#include "crosslink_species.hpp"

void CrosslinkSpecies::Init(system_parameters *params,
                            species_base_parameters *sparams,
                            space_struct *space) {
  Species::Init(params, sparams, space);
  k_on_ = sparams_.k_on;
  k_off_ = sparams_.k_off;
  xlink_concentration_ = sparams_.concentration;
  sparams_.num = (int) floor(sparams_.concentration * space->volume);
}

void CrosslinkSpecies::AddMember() {
  Species::AddMember();
  members_.back().InitInteractionEnvironment(&lut_);
}

void CrosslinkSpecies::InitInteractionEnvironment(std::vector<Object *> *objs,
                                                  double *obj_vol,
                                                  bool *update) {
  objs_ = objs;
  obj_volume_ = obj_vol;
  update_ = update;
  /* TODO Lookup table only works for filament objects. Generalize? */
  lut_.Init(sparams_.k_spring / 2, sparams_.rest_length, 1);
}

void CrosslinkSpecies::CalculateBindingFree() {
  /* Check crosslink binding */
  double free_concentration = (sparams_.num - n_members_) / space_->volume;
  if (gsl_rng_uniform_pos(rng_.r) <=
      free_concentration * (*obj_volume_) * k_on_ * params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
    *update_ = true;
  }
}

/* Returns a random object with selection probability proportional to object
   volume */
Object *CrosslinkSpecies::GetRandomObject() {
  double roll = (*obj_volume_) * gsl_rng_uniform_pos(rng_.r);
  double vol = 0;
  for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
    vol += (*obj)->GetVolume();
    if (vol > roll) {
#ifdef TRACE
      Logger::Trace("Binding free crosslink to random object: xl %d -> obj %d",
                    members_.back().GetOID(), (*obj)->GetOID());
#endif
      return *obj;
    }
  }
  Logger::Error("CrosslinkSpecies::GetRandomObject should never get here!");
}

/* A crosslink binds to an object from solution */
void CrosslinkSpecies::BindCrosslink() {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  AddMember();
  members_.back().AttachObjRandom(GetRandomObject());
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to bonds upon a checkpoint reload */
void CrosslinkSpecies::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetAnchors(ixors);
  }
}

void CrosslinkSpecies::UpdatePositions() {
  /* Only do this every other step (assuming flexible filaments with midstep)
   */
  if (params_->i_step % 2 == 0) {
    /* First update bound crosslinks state and positions */
    UpdateBoundCrosslinks();
    /* Calculate implicit binding of crosslinks from solution */
    CalculateBindingFree();
  } else {
    /* Apply tether forces from doubly-bound crosslinks onto anchored objects.
       We do this every half step only, because the fullstep tether forces are
       handled by the crosslink update on even steps */
    ApplyCrosslinkTetherForces();
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinks() {
  n_members_ = 0;
  /* Update anchor positions to their attached meshes and calculate anchor
     forces */
  UpdateBoundCrosslinkForces();
  /* Apply anchor forces on bound objects sequentially */
  ApplyCrosslinkTetherForces();
  /* Update anchor positions from diffusion, walking */
  UpdateBoundCrosslinkPositions();
  /* Remove crosslinks that came unbound */
  members_.erase(std::remove_if(members_.begin(), members_.end(),
                                [](Crosslink x) { return x.IsUnbound(); }),
                 members_.end());
  /* Get the number of bound crosslinks so we know what the current
     concentration of free crosslinks is */
  n_members_ = members_.size();
}

/* This must be done sequentially to avoid racy conditions when accessing
   bound object's forces */
void CrosslinkSpecies::ApplyCrosslinkTetherForces() {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->ApplyTetherForces();
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinkForces() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        xlink->UpdateCrosslinkForces();
        if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    xlink->UpdateCrosslinkForces();
    if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}
void CrosslinkSpecies::UpdateBoundCrosslinkPositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        xlink->UpdateCrosslinkPositions();
        /* Xlink is no longer bound, return to solution */
        if (xlink->IsUnbound()) {
          *update_ = true;
          /* If a crosslink enters or leaves the singly state, we need to
           * update xlink interactors */
        } else if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    xlink->UpdateCrosslinkPositions();
    /* Xlink is no longer bound, return to solution */
    if (xlink->IsUnbound()) {
      *update_ = true;
      /* If a crosslink enters or leaves the singly state, we need to update
       * xlink interactors */
    } else if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}

void CrosslinkSpecies::CleanUp() { members_.clear(); }

void CrosslinkSpecies::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkSpecies::ReadSpecs() {
  if (ispec_file_.eof()) {
    Logger::Info("EOF reached for spec file in CrosslinkSpecies");
    early_exit = true;
    return;
  }
  if (!ispec_file_.is_open()) {
    Logger::Warning("ERROR. Spec file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  n_members_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_members_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_members_ == -1) {
    Logger::Info("EOF reached for spec file in CrosslinkSpecies");
    early_exit = true;
    return;
  }
  if (n_members_ == 0) {
    members_.clear();
  } else if (n_members_ != members_.size()) {
    Crosslink xlink;
    xlink.Init(&sparams_);
    xlink.InitInteractionEnvironment(&lut_);
    members_.resize(n_members_, xlink);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadSpec(ispec_file_);
  }
}
