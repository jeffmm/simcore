#include "crosslink_species.hpp"

CrosslinkSpecies::CrosslinkSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::crosslink);
}
void CrosslinkSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  k_on_ = sparams_.k_on;
  k_off_ = sparams_.k_off;
  xlink_concentration_ = sparams_.concentration;
  sparams_.num = (int)round(sparams_.concentration * space_->volume);
}

void CrosslinkSpecies::AddMember() {
  Species::AddMember();
  members_.back().InitInteractionEnvironment(&lut_);
  *update_ = true;
}

void CrosslinkSpecies::InitInteractionEnvironment(std::vector<Object *> *objs,
                                                  double *obj_vol,
                                                  bool *update) {
  objs_ = objs;
  obj_volume_ = obj_vol;
  update_ = update;
  /* TODO Lookup table only works for filament objects. Generalize? */
  lut_.Init(sparams_.k_spring / 2, sparams_.rest_length, 1);
  /* Integral cutoff */
  double small = 1e-4; // Value defined in lookup_table
  sparams_.r_capture =
      sqrt(-2 * log(small) / sparams_.k_spring) + 1 + sparams_.rest_length;
}

void CrosslinkSpecies::InsertCrosslinks() {
  // Random insertion (default) implies crosslinks in solution
  if (sparams_.insertion_type.compare("random") == 0) {
    sparams_.static_flag = false;
    return;
  } else if (sparams_.insertion_type.compare("centered") == 0) {
    sparams_.num = 1;
    sparams_.static_flag = true;
    AddMember();
    double pos[3] = {0, 0, 0};
    double u[3] = {0, 0, 0};
    u[params_->n_dim - 1] = 1.0;
    members_.back().InsertAt(pos, u);
  } else if (sparams_.insertion_type.compare("random_grid") == 0) {
    sparams_.static_flag = true;
    if (params_->n_dim == 3) {
      if (space_->type == +boundary_type::none ||
          space_->type == +boundary_type::box) {
        sparams_.num = (int)round(4 * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::sphere) {
        sparams_.num = (int)round(M_PI * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::budding) {
        double R = space_->radius;
        double r = space_->bud_radius;
        double d = space_->bud_height;
        sparams_.num = (int)round(
            M_PI * SQR(R) + M_PI * SQR(r) -
            SQR(r) * acos((SQR(d) + SQR(r) - SQR(R)) / (2 * d * r)) -
            SQR(R) * acos((SQR(d) + SQR(R) - SQR(r)) / (2 * d * R)) +
            0.5 * sqrt((R + r - d) * (r + d - R) * (R + d - r) * (R + r + d)));
      } else {
        Logger::Error("Boundary type not recognized in CrosslinkSpecies");
      }
    }
    for (int i = 0; i < sparams_.num; ++i) {
      AddMember();
      members_.back().InsertRandom();
      /* If in 3D, zero out the third dimension so the anchor is on a plane */
      if (params_->n_dim == 3) {
        double projected_pos[3] = {0};
        double same_u[3] = {0};
        const double *const pos = members_.back().GetPosition();
        const double *const u = members_.back().GetOrientation();
        for (int j = 0; j < 3; ++j) {
          projected_pos[j] = pos[j];
          same_u[j] = u[j];
        }
        projected_pos[2] = 0;
        members_.back().InsertAt(projected_pos, same_u);
      }
    }
  } else if (sparams_.insertion_type.compare("random_boundary") == 0) {
    if (space_->type == +boundary_type::none) {
      Logger::Error("Crosslinker insertion type \"random boundary\" requires a"
                    " boundary for species insertion.");
    } else if (space_->type == +boundary_type::sphere) {
      if (params_->n_dim == 2) {
        sparams_.num =
            (int)round(2.0 * M_PI * space_->radius * xlink_concentration_);
      } else {
        sparams_.num =
            (int)round(4 * M_PI * SQR(space_->radius) * xlink_concentration_);
      }
      double pos[3] = {0};
      double u[3] = {0};
      for (int i = 0; i < sparams_.num; ++i) {
        rng_.RandomUnitVector(params_->n_dim, u);
        for (int j = 0; j < params_->n_dim; ++j) {
          pos[j] = u[j] * space_->radius;
        }
        AddMember();
        members_.back().InsertAt(pos, u);
      }
    } else if (space_->type == +boundary_type::budding){

    } else {
      Logger::Error("Boundary type not recognized in CrosslinkSpecies");
    }
  } else {
    Logger::Error("Insertion type %s not implemented yet for crosslinks",
                  sparams_.insertion_type.c_str());
  }
}

void CrosslinkSpecies::CalculateBindingFree() {
  /* Static crosslinks are never free */
  if (sparams_.static_flag) {
    return;
  }
  /* Check crosslink binding */
  double free_concentration = (sparams_.num - n_members_) / space_->volume;
  if (rng_.RandomUniform() <=
      free_concentration * (*obj_volume_) * k_on_ * params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
  }
}

/* Returns a random object with selection probability proportional to object
   volume */
Object *CrosslinkSpecies::GetRandomObject() {
  double roll = (*obj_volume_) * rng_.RandomUniform();
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
  // for (auto it=members_.begin(); it!=members_.end(); ++it) {
  // it->SanityCheck();
  //}
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
  if (!sparams_.static_flag) {
    members_.erase(std::remove_if(members_.begin(), members_.end(),
                                  [](Crosslink x) { return x.IsUnbound(); }),
                   members_.end());
  }
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
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
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
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
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
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
        xlink->UpdateCrosslinkPositions();
        /* Xlink is no longer bound, return to solution */
        if (xlink->IsUnbound()) {
          if (sparams_.static_flag) {
            Logger::Error("Static crosslinks became unbound");
          }
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
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
    xlink->UpdateCrosslinkPositions();
    /* Xlink is no longer bound, return to solution */
    if (xlink->IsUnbound()) {
      if (sparams_.static_flag) {
        Logger::Error("Static crosslinks became unbound");
      }
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
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
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
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (n_members_ == 0) {
    members_.clear();
  } else if (n_members_ != members_.size()) {
    Crosslink xlink(rng_.GetSeed());
    xlink.Init(&sparams_);
    xlink.InitInteractionEnvironment(&lut_);
    xlink.SetSID(GetSID());
    members_.resize(n_members_, xlink);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadSpec(ispec_file_);
  }
}

const double CrosslinkSpecies::GetConcentration() const {
  return sparams_.concentration;
}
const double CrosslinkSpecies::GetRCutoff() const { return sparams_.r_capture; }
