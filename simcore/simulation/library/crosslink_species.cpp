#include "crosslink_species.hpp"

void CrosslinkSpecies::Init(system_parameters *params,
                            species_base_parameters *sparams,
                            space_struct *space) {
  Species::Init(params, sparams, space);
  k_on_ = sparams_.k_on;
  k_off_ = sparams_.k_off;
  xlink_concentration_ = sparams_.concentration;
}

void CrosslinkSpecies::AddMember() {
  Crosslink newmember;
  Logger::Trace("Adding member to species %s, member number %d, member id %d",
      GetSID()._to_string(), n_members_+1, newmember.GetOID());
  members_.push_back(newmember);
  members_.back().SetSID(GetSID());
  members_.back().Init(&sparams_);
  members_.back().InitInteractionEnvironment(&lut_);
  /* Keep track of bound bound anchors, bound crosslinks, and
   * concentration of free crosslinks */
  n_members_++;
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
  double concentration = xlink_concentration_ - n_members_ / space_->volume;
  if (gsl_rng_uniform_pos(rng_.r) <=
      concentration * (*obj_volume_) * k_on_ * params_->delta) {
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

void CrosslinkSpecies::UpdateCrosslinks() {
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

void CrosslinkSpecies::WriteSpecs() {
  /* Write the vector sizes, singly then doubly */
  n_members_ = members_.size();
  ospec_file_.write(reinterpret_cast<char *>(&n_members_), sizeof(int));

  /* Write individual crosslink specs, first singly then doubly */
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->WriteSpec(ospec_file_);
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

void CrosslinkSpecies::WriteCheckpoints() {
  /* Try to open the file */
  std::fstream ocheck_file(checkpoint_file_, std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    Logger::Error("Output file %s did not open\n", checkpoint_file_.c_str());
  }

  /* Write RNG state */
  long seed = rng_.GetSeed();
  ocheck_file.write(reinterpret_cast<char *>(&seed), sizeof(long));
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck_file.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck_file.write(reinterpret_cast<char *>(rng_state), rng_size);

  /* Write crosslink vector sizes, singly then doubly */
  n_members_ = members_.size();
  ocheck_file.write(reinterpret_cast<char *>(&n_members_), sizeof(int));

  /* Write crosslink checkpoints, singly then doubly */
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->WriteCheckpoint(ocheck_file);
  }

  /* Close the file */
  ocheck_file.close();
}

void CrosslinkSpecies::ReadCheckpoints() {
  /* Try to open the file */
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    Logger::Error("Output file %s did not open\n", checkpoint_file_.c_str());
  }

  /* Read RNG state */
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  long seed = -1;
  icheck_file.read(reinterpret_cast<char *>(&seed), sizeof(long));
  icheck_file.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck_file.read(reinterpret_cast<char *>(rng_state), rng_size);
  /* Read xlink vector sizes, first singly then doubly */
  n_members_ = -1;
  icheck_file.read(reinterpret_cast<char *>(&n_members_), sizeof(int));

  /* Prepare the xlink vectors */
  Crosslink xlink;
  xlink.Init(&sparams_);
  xlink.InitInteractionEnvironment(&lut_);
  members_.resize(n_members_, xlink);

  /* Read the crosslink checkpoints */
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadCheckpoint(icheck_file);
  }
  /* Close the file */
  icheck_file.close();
  rng_.SetSeed(seed);
}

void CrosslinkSpecies::InitOutputFiles() {
  if (GetSpecFlag())
    InitSpecFile();
  if (GetCheckpointFlag())
    InitCheckpoints();
}

void CrosslinkSpecies::InitSpecFile() {
  std::string spec_file_name =
      params_->run_name + "_xlink_" + GetSpeciesName() + ".spec";
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary);
  if (!ospec_file_.is_open()) {
    Logger::Error("Output file %s did not open\n", spec_file_name.c_str());
  }
  int n_spec = GetNSpec();
  ospec_file_.write(reinterpret_cast<char *>(&params_->n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&n_spec), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&params_->delta), sizeof(double));
}

void CrosslinkSpecies::InitCheckpoints() {
  checkpoint_file_ =
      params_->run_name + "_xlink_" + GetSpeciesName() + ".checkpoint";
}

bool CrosslinkSpecies::InitSpecFileInputFromFile(std::string spec_file_name) {
  ispec_file_.open(spec_file_name, std::ios::in | std::ios::binary);
  if (!ispec_file_.is_open()) {
    return false;
  }
  // long n_steps;
  int n_spec, n_steps;
  double delta;
  ispec_file_.read(reinterpret_cast<char *>(&n_steps), sizeof(int));
  ispec_file_.read(reinterpret_cast<char *>(&n_spec), sizeof(int));
  ispec_file_.read(reinterpret_cast<char *>(&delta), sizeof(double));
  if (n_steps != params_->n_steps || n_spec != GetNSpec() ||
      delta != params_->delta) {
    Logger::Warning("Input file %s does not match parameter file\n",
                    "n_steps: %d %d, n_spec: %d %d, delta: %2.2f %2.2f",
                    spec_file_name.c_str(), n_steps, params_->n_steps, n_spec,
                    GetNSpec(), delta, params_->delta);
  }
  ReadSpecs();
  return true;
}

void CrosslinkSpecies::InitSpecFileInput() {
  std::string spec_file_name =
      params_->run_name + "_xlink_" + GetSpeciesName() + ".spec";
  if (!InitSpecFileInputFromFile(spec_file_name)) {
    Logger::Error("Input file %s did not open", spec_file_name.c_str());
  }
}

// void CrosslinkSpecies::UpdatePeriodic() {
// for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
// xlink->UpdatePeriodic();
//}
//}

void CrosslinkSpecies::LoadFromCheckpoints() {
  checkpoint_file_ =
      params_->checkpoint_run_name + "_xlink_" + GetSpeciesName() + ".checkpoint";
  if (!GetCheckpointFlag()) {
    Logger::Error("Checkpoint file %s not available for parameter file!",
                  checkpoint_file_.c_str());
  }
  ReadCheckpoints();
  InitOutputFiles();
}

void CrosslinkSpecies::InitOutputs(bool reading_inputs, bool reduce_flag,
                                   bool with_reloads) {
  if (sparams_.concentration < 1e-12) {
    return;
  }
  if (params_->load_checkpoint) {
    LoadFromCheckpoints();
  } else if (reading_inputs) {
    if (with_reloads) {
      Logger::Warning("Crosslinks do not yet support reload functionality.");
    }
    InitSpecFileInput();
    if (reduce_flag) {
      Logger::Warning("Crosslinks do not yet support reduce functionality.");
    } else if (params_->checkpoint_from_spec) {
      InitCheckpoints();
    }
  } else {
    InitOutputFiles();
  }
}
void CrosslinkSpecies::WriteOutputs() {
  if (sparams_.concentration < 1e-12) {
    return;
  }
  if (GetSpecFlag() && (params_->i_step % GetNSpec() == 0)) {
    WriteSpecs();
  }
  if (GetCheckpointFlag() && (params_->i_step % GetNCheckpoint() == 0)) {
    WriteCheckpoints();
  }
}

void CrosslinkSpecies::ReadInputs() {
  if (sparams_.concentration < 1e-12) {
    return;
  }
  if (GetSpecFlag() && (params_->i_step % GetNSpec() == 0)) {
    ReadSpecs();
    if (params_->checkpoint_from_spec) {
      WriteCheckpoints();
    }
  }
}
