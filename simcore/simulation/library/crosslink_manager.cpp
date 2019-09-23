#include "crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, space_struct *space,
                            MinimumDistance *mindist,
                            std::vector<Object *> *objs) {
  params_ = params;
  space_ = space;
  mindist_ = mindist;
  objs_ = objs;
  k_on_ = params_->crosslink.k_on;
  k_off_ = params_->crosslink.k_off;
  xlink_concentration_ = params_->crosslink.concentration;
  obj_volume_ = 0;
  n_xlinks_ = 0;
  n_anchors_bound_ = 0;
  n_spec_ = params_->crosslink.n_spec;
  n_checkpoint_ = params_->crosslink.n_checkpoint;
  checkpoint_flag_ = params_->crosslink.checkpoint_flag;
  spec_flag_ = params_->crosslink.spec_flag;
  update_ = false;
  /* TODO Lookup table only works for filament objects. Generalize? */
  lut_.Init(params_->crosslink.k_spring / 2, params_->crosslink.rest_length,
            params_->filament.diameter);
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsVolume() {
  obj_volume_ = 0;
  for (auto it = objs_->begin(); it != objs_->end(); ++it) {
    obj_volume_ += (*it)->GetVolume();
  }
}

/* Whether to reinsert anchors into the interactors list */
bool CrosslinkManager::CheckUpdate() {
  if (update_) {
    update_ = false;
    return true;
  }
  return false;
}

void CrosslinkManager::CalculateBindingFree() {
  /* Check crosslink binding */
  UpdateObjsVolume();
  double concentration = xlink_concentration_ - n_xlinks_ / space_->volume;
  if (gsl_rng_uniform_pos(rng_.r) <=
      concentration * obj_volume_ * k_on_ * params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
    update_ = true;
  }
}

/* Returns a random object with selection probability proportional to object
   volume */
Object *CrosslinkManager::GetRandomObject() {
  double roll = obj_volume_ * gsl_rng_uniform_pos(rng_.r);
  double vol = 0;
  for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
    vol += (*obj)->GetVolume();
    if (vol > roll)
      return *obj;
  }
  error_exit("CrosslinkManager::GetRandomObject should never get here!");
}

/* A crosslink binds to an object from solution */
void CrosslinkManager::BindCrosslink() {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  Crosslink xl;
  xlinks_.push_back(xl);
  xlinks_.back().Init(mindist_, &lut_);
  xlinks_.back().AttachObjRandom(GetRandomObject());
  /* Keep track of bound bound anchors, bound crosslinks, and
   * concentration of free crosslinks */
  n_anchors_bound_++;
  n_xlinks_++;
}

void CrosslinkManager::GetInteractors(std::vector<Object *> *ixors) {
  std::vector<Object *> ix;
  /* Clear the crosslink neighbors as we go, since we need to repopulate the
     neighbor lists if we are repopulating the interactor vector */
  for (auto xlink = xlinks_.begin(); xlink != xlinks_.end(); ++xlink) {
    xlink->GetInteractors(&ix);
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void CrosslinkManager::GetAnchorInteractors(std::vector<Object *> *ixors) {
  std::vector<Object *> ix;
  for (auto xlink = xlinks_.begin(); xlink != xlinks_.end(); ++xlink) {
    xlink->GetAnchors(&ix);
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void CrosslinkManager::UpdateCrosslinks() {
  update_ = false;
  /* Only do this every other step (assuming flexible filaments with midstep) */
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

void CrosslinkManager::ApplyCrosslinkTetherForces() {
  for (auto xlink = xlinks_.begin(); xlink != xlinks_.end(); ++xlink) {
    xlink->ApplyTetherForces();
  }
}

void CrosslinkManager::UpdateBoundCrosslinks() {
  n_anchors_bound_ = 0;
  n_xlinks_ = 0;
  for (auto xlink = xlinks_.begin(); xlink != xlinks_.end(); ++xlink) {
    bool init_state = xlink->IsSingly();
    n_xlinks_++;
    xlink->UpdateCrosslink();
    /* Xlink is no longer bound, return to solution */
    if (xlink->IsUnbound()) {
      update_ = true;
      n_xlinks_--;
    } else if (xlink->IsSingly()) {
      n_anchors_bound_++;
    } else if (xlink->IsDoubly()) {
      n_anchors_bound_ += 2;
    }
    /* If a crosslink enters or leaves the singly state, we need to update
     * xlink interactors */
    if (xlink->IsSingly() != init_state) {
      update_ = true;
    }
  }

  /* Remove unbound crosslinks */
  xlinks_.erase(std::remove_if(xlinks_.begin(), xlinks_.end(),
                               [](Crosslink x) { return x.IsUnbound(); }),
                xlinks_.end());
}

void CrosslinkManager::Clear() { xlinks_.clear(); }

void CrosslinkManager::Draw(std::vector<graph_struct *> *graph_array) {
  for (auto it = xlinks_.begin(); it != xlinks_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToAnchor(Object *anchor, Object *neighbor) {
  if (anchor->GetSID() != +species_id::crosslink) {
    error_exit(
        "BindCrosslinkObj expected crosslink object, got generic object.");
  }
  Anchor *a = dynamic_cast<Anchor *>(anchor);
  a->AddNeighbor(neighbor);
}

void CrosslinkManager::WriteSpecs() {
  /* Write the vector sizes, singly then doubly */
  n_xlinks_ = xlinks_.size();
  ospec_file_.write(reinterpret_cast<char *>(&n_xlinks_), sizeof(int));

  /* Write individual crosslink specs, first singly then doubly */
  for (auto it = xlinks_.begin(); it != xlinks_.end(); ++it) {
    it->WriteSpec(ospec_file_);
  }
}

void CrosslinkManager::ReadSpecs() {
  if (ispec_file_.eof()) {
    printf("  EOF reached\n");
    early_exit = true;
    return;
  }
  if (!ispec_file_.is_open()) {
    printf(" ERROR: Spec file unexpectedly not open! Exiting early.\n");
    early_exit = true;
    return;
  }
  n_xlinks_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_xlinks_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_xlinks_ == -1) {
    printf("  EOF reached in species\n");
    early_exit = true;
    return;
  }
  if (n_xlinks_ == 0) {
    xlinks_.clear();
  } else if (n_xlinks_ != xlinks_.size()) {
    Crosslink xlink;
    xlink.Init(mindist_, &lut_);
    xlinks_.resize(n_xlinks_, xlink);
  }
  n_anchors_bound_ = 0;
  for (auto it = xlinks_.begin(); it != xlinks_.end(); ++it) {
    it->ReadSpec(ispec_file_);
    n_anchors_bound_++;
    if (it->IsDoubly()) {
      n_anchors_bound_++;
    }
  }
}

void CrosslinkManager::WriteCheckpoints() {
  /* Try to open the file */
  std::fstream ocheck_file(checkpoint_file_, std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    error_exit("Output file %s did not open\n", checkpoint_file_.c_str());
  }

  /* Write RNG state */
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck_file.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck_file.write(reinterpret_cast<char *>(rng_state), rng_size);

  /* Write crosslink vector sizes, singly then doubly */
  n_xlinks_ = xlinks_.size();
  ocheck_file.write(reinterpret_cast<char *>(&n_xlinks_), sizeof(int));

  /* Write crosslink checkpoints, singly then doubly */
  for (auto it = xlinks_.begin(); it != xlinks_.end(); ++it) {
    it->WriteCheckpoint(ocheck_file);
  }

  /* Close the file */
  ocheck_file.close();
}

void CrosslinkManager::ReadCheckpoints() {
  /* Try to open the file */
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    error_exit("Output file %s did not open\n", checkpoint_file_.c_str());
  }

  /* Read RNG state */
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck_file.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck_file.read(reinterpret_cast<char *>(rng_state), rng_size);

  /* Read xlink vector sizes, first singly then doubly */
  n_xlinks_ = -1;
  icheck_file.read(reinterpret_cast<char *>(&n_xlinks_), sizeof(int));

  /* Prepare the xlink vectors */
  Crosslink xlink;
  xlink.Init(mindist_, &lut_);
  xlinks_.resize(n_xlinks_, xlink);

  /* Read the crosslink checkpoints */
  n_anchors_bound_ = 0;
  for (auto it = xlinks_.begin(); it != xlinks_.end(); ++it) {
    it->ReadCheckpoint(icheck_file);
    n_anchors_bound_++;
    if (it->IsDoubly()) {
      n_anchors_bound_++;
    }
  }
  /* Close the file */
  icheck_file.close();
}

void CrosslinkManager::InitOutputFiles() {
  if (params_->crosslink.spec_flag)
    InitSpecFile();
  if (params_->crosslink.checkpoint_flag)
    InitCheckpoints();
}

void CrosslinkManager::InitSpecFile() {
  std::string spec_file_name = params_->run_name + "_crosslink.spec";
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary);
  if (!ospec_file_.is_open()) {
    error_exit("Output file %s did not open\n", spec_file_name.c_str());
  }
  ospec_file_.write(reinterpret_cast<char *>(&params_->n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&params_->crosslink.n_spec),
                    sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&params_->delta), sizeof(double));
}

void CrosslinkManager::InitCheckpoints() {
  checkpoint_file_ = params_->run_name + "_crosslink.checkpoint";
}

bool CrosslinkManager::InitSpecFileInputFromFile(std::string spec_file_name) {
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
  if (n_steps != params_->n_steps || n_spec != params_->crosslink.n_spec ||
      delta != params_->delta) {
    printf("\nn_steps: %d %d, ", n_steps, params_->n_steps);
    printf("n_spec: %d %d, ", n_spec, params_->crosslink.n_spec);
    printf("delta: %2.2f %2.2f\n", delta, params_->delta);
    warning("Input file %s does not match parameter file\n",
            spec_file_name.c_str());
  }
  ReadSpecs();
  return true;
}

void CrosslinkManager::InitSpecFileInput() {
  std::string spec_file_name = params_->run_name + "_crosslink.spec";
  if (!InitSpecFileInputFromFile(spec_file_name)) {
    error_exit("Input file %s did not open", spec_file_name.c_str());
  }
}

void CrosslinkManager::LoadFromCheckpoints() {
  checkpoint_file_ = params_->checkpoint_run_name + "_crosslink.checkpoint";
  if (!params_->crosslink.checkpoint_flag) {
    error_exit("Checkpoint file %s not available for parameter file!",
               checkpoint_file_.c_str());
  }
  ReadCheckpoints();
  InitOutputFiles();
}

void CrosslinkManager::InitOutputs(bool reading_inputs, bool reduce_flag,
                                   bool with_reloads) {
  if (params_->load_checkpoint) {
    LoadFromCheckpoints();
  } else if (reading_inputs) {
    if (with_reloads) {
      warning("Crosslinks do not yet support reload functionality.");
    }
    InitSpecFileInput();
    if (reduce_flag) {
      warning("Crosslinks do not yet support reduce functionality.");
    } else if (params_->checkpoint_from_spec) {
      InitCheckpoints();
    }
  } else {
    InitOutputFiles();
  }
}
void CrosslinkManager::WriteOutputs() {
  if (spec_flag_ && (params_->i_step % n_spec_ == 0)) {
    WriteSpecs();
  }
  if (checkpoint_flag_ && (params_->i_step % n_checkpoint_ == 0)) {
    WriteCheckpoints();
  }
}

void CrosslinkManager::ReadInputs() {
  if (spec_flag_ && (params_->i_step % n_spec_ == 0)) {
    ReadSpecs();
    if (params_->checkpoint_from_spec) {
      WriteCheckpoints();
    }
  }
}
