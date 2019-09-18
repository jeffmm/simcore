#include "crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, space_struct *space,
                            MinimumDistance *mindist,
                            std::vector<Object *> *objs) {
  params_ = params;
  space_ = space;
  mindist_ = mindist;
  objs_ = objs;
  rng_.Init();
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
  /* TODO, keep track of individual object binding probabilities,
   * which would end up being individual object volume/total
   * object volume. For now, will work correctly with equally-sized
   * objects. */
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
  // Check crosslink unbinding
  // if (gsl_rng_uniform_pos(rng_.r) <= k_off_*n_anchors_bound_*params_->delta)
  // {
  //[> Remove a random anchor from the system <]
  // UnbindCrosslink();
  // update_ = true;
  //}
  /* Check crosslink binding */
  double concentration = xlink_concentration_ - n_xlinks_ / space_->volume;
  if (gsl_rng_uniform_pos(rng_.r) <=
      concentration * obj_volume_ * k_on_ * params_->delta) {
    /* Create a new crosslink and bind an anchor to a random object
     * in the system */
    BindCrosslink();
    update_ = true;
  }
}

/* A crosslink binds to an object from solution */
void CrosslinkManager::BindCrosslink() {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  Crosslink xl;
  xlinks_singly_.push_back(xl);
  xlinks_singly_.back().Init(mindist_, &lut_);
  /* Attach to random object in system */
  /* TODO Should weight probability of selecting object
     by object volume in the case of different sized
     objects */
  int i_obj = gsl_rng_uniform_int(rng_.r, objs_->size());
  xlinks_singly_.back().AttachObjRandom((*objs_)[i_obj]);
  /* Keep track of bound bound anchors, bound crosslinks, and
   * concentration of free crosslinks */
  n_anchors_bound_++;
  n_xlinks_++;
  // printf("%d \n", n_anchors_bound_);
}

/* An unbinding event for a single, random anchor in the system */
void CrosslinkManager::UnbindCrosslink() {
  int i_anchor = gsl_rng_uniform_int(rng_.r, n_anchors_bound_);
  /* If i_anchor < the number of singly-bound anchors, unbind a
   * singly-bound anchor */
  int n_singly = xlinks_singly_.size();
  if (i_anchor < n_singly) {
    RemoveCrosslink(i_anchor);
  } else {
    i_anchor -= n_singly;
    int i_doubly = i_anchor / 2;
    int i_which = i_anchor % 2;
    xlinks_doubly_[i_doubly].UnbindAnchor(i_which);
    DoublyToSingly(i_doubly);
  }
  n_anchors_bound_--;
}

void CrosslinkManager::RemoveCrosslink(int i_xlink) {
  if (xlinks_singly_.empty()) {
    return;
  }
  if (xlinks_singly_.size() > 1) {
    auto it = (xlinks_singly_.begin() + i_xlink);
    it->UnbindAnchor();
    std::iter_swap(xlinks_singly_.begin() + i_xlink, xlinks_singly_.end() - 1);
    xlinks_singly_.pop_back();
  } else {
    xlinks_singly_.clear();
  }
  n_xlinks_--;
  update_ = true;
}

void CrosslinkManager::DoublyToSingly(int i_doubly) {
  if (xlinks_doubly_.size() > 1) {
    std::iter_swap(xlinks_doubly_.begin() + i_doubly, xlinks_doubly_.end() - 1);
    xlinks_singly_.push_back(*(xlinks_doubly_.end() - 1));
    xlinks_doubly_.pop_back();
  } else {
    xlinks_singly_.push_back(*(xlinks_doubly_.begin()));
    xlinks_doubly_.clear();
  }
  xlinks_singly_.back().SetSingly();
}

void CrosslinkManager::GetInteractors(std::vector<Object *> *ixors) {
  ixors->clear();
  std::vector<Object *> ix;
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end();
       ++xlink) {
    ix.push_back(&(*xlink));
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void CrosslinkManager::GetAnchorInteractors(std::vector<Object *> *ixors) {
  ixors->clear();
  std::vector<Object *> ix;
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end();
       ++xlink) {
    xlink->GetAnchors(&ix);
  }
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end();
       ++xlink) {
    xlink->GetAnchors(&ix);
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void CrosslinkManager::UpdateCrosslinks() {
  /* First update bound crosslinks, then update binding */
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end();
       ++xlink) {
    xlink->UpdateCrosslink();
  }
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end();
       ++xlink) {
    xlink->UpdateCrosslink();
  }
  /* Check if we had any crosslinking events */
  for (auto xlink = xlinks_singly_.begin(); xlink != xlinks_singly_.end();
       ++xlink) {
    if (xlink->IsDoubly()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_singly_.end() - 1);
      xlinks_doubly_.push_back(*(xlinks_singly_.end() - 1));
      xlinks_singly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_++;
    } else if (xlink->IsUnbound()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_singly_.end() - 1);
      xlinks_singly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_--;
      n_xlinks_--;
    }
  }
  /* Check if we had any crosslinks break */
  for (auto xlink = xlinks_doubly_.begin(); xlink != xlinks_doubly_.end();
       ++xlink) {
    if (xlink->IsSingly()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_doubly_.end() - 1);
      xlinks_singly_.push_back(*(xlinks_doubly_.end() - 1));
      xlinks_doubly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_--;
    } else if (xlink->IsUnbound()) {
      auto it = xlink;
      std::iter_swap(it, xlinks_singly_.end() - 1);
      xlinks_singly_.pop_back();
      xlink--;
      update_ = true;
      n_anchors_bound_ -= 2;
      n_xlinks_--;
    }
  }

  /* Calculate binding of free crosslinks (in solution)
     and unbinding of bound crosslinks */
  CalculateBindingFree();
}

void CrosslinkManager::Clear() {
  xlinks_singly_.clear();
  xlinks_doubly_.clear();
}

void CrosslinkManager::Draw(std::vector<graph_struct *> *graph_array) {
  for (auto it = xlinks_singly_.begin(); it != xlinks_singly_.end(); ++it) {
    it->Draw(graph_array);
  }
  for (auto it = xlinks_doubly_.begin(); it != xlinks_doubly_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToXlink(Object *xlink, Object *neighbor) {
  if (xlink->GetSID() != +species_id::crosslink) {
    error_exit(
        "BindCrosslinkObj expected crosslink object, got generic object.");
  }
  Crosslink *xl = dynamic_cast<Crosslink *>(xlink);
  xl->AddNeighbor(neighbor);
}

void CrosslinkManager::WriteSpecs() {
  /* Write the vector sizes, singly then doubly */
  int size_singly = xlinks_singly_.size();
  int size_doubly = xlinks_doubly_.size();
  ospec_file_.write(reinterpret_cast<char *>(&size_singly), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&size_doubly), sizeof(int));

  /* Write individual crosslink specs, first singly then doubly */
  for (auto it = xlinks_singly_.begin(); it != xlinks_singly_.end(); ++it) {
    it->WriteSpec(ospec_file_);
  }
  for (auto it = xlinks_doubly_.begin(); it != xlinks_doubly_.end(); ++it) {
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
  int size_singly = -1;
  int size_doubly = -1;
  ispec_file_.read(reinterpret_cast<char *>(&size_singly), sizeof(int));
  ispec_file_.read(reinterpret_cast<char *>(&size_doubly), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (size_singly == -1 || size_doubly == -1) {
    printf("  EOF reached in species\n");
    early_exit = true;
    return;
  }
  if (size_singly == 0) {
    xlinks_singly_.clear();
  } else if (size_singly != xlinks_singly_.size()) {
    Crosslink xlink;
    xlink.Init(mindist_, &lut_);
    xlinks_singly_.resize(size_singly, xlink);
  }
  if (size_doubly == 0) {
    xlinks_doubly_.clear();
  } else if (size_doubly != xlinks_doubly_.size()) {
    Crosslink xlink;
    xlink.Init(mindist_, &lut_);
    xlink.SetDoubly();
    xlinks_doubly_.resize(size_doubly, xlink);
  }
  n_xlinks_ = 0;
  n_anchors_bound_ = 0;
  for (auto it = xlinks_singly_.begin(); it != xlinks_singly_.end(); ++it) {
    // it->Init(mindist_, &lut_);
    it->ReadSpec(ispec_file_);
    n_xlinks_++;
    n_anchors_bound_++;
  }
  for (auto it = xlinks_doubly_.begin(); it != xlinks_doubly_.end(); ++it) {
    // it->Init(mindist_, &lut_);
    it->ReadSpec(ispec_file_);
    n_xlinks_++;
    n_anchors_bound_ += 2;
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
  int size_singly = xlinks_singly_.size();
  int size_doubly = xlinks_doubly_.size();
  ocheck_file.write(reinterpret_cast<char *>(&size_singly), sizeof(int));
  ocheck_file.write(reinterpret_cast<char *>(&size_doubly), sizeof(int));

  /* Write crosslink checkpoints, singly then doubly */
  for (auto it = xlinks_singly_.begin(); it != xlinks_singly_.end(); ++it) {
    it->WriteCheckpoint(ocheck_file);
  }
  for (auto it = xlinks_doubly_.begin(); it != xlinks_doubly_.end(); ++it) {
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
  int size_singly = 0;
  int size_doubly = 0;
  icheck_file.read(reinterpret_cast<char *>(&size_singly), sizeof(int));
  icheck_file.read(reinterpret_cast<char *>(&size_doubly), sizeof(int));

  /* Prepare the xlink vectors */
  Crosslink xlink;
  xlink.Init(mindist_, &lut_);
  xlinks_singly_.resize(size_singly, xlink);
  xlinks_doubly_.resize(size_doubly, xlink);

  /* Read the crosslink checkpoints */
  for (auto it = xlinks_singly_.begin(); it != xlinks_singly_.end(); ++it) {
    // it->Init(mindist_, &lut_);
    it->ReadCheckpoint(icheck_file);
    n_xlinks_++;
    n_anchors_bound_++;
  }
  for (auto it = xlinks_doubly_.begin(); it != xlinks_doubly_.end(); ++it) {
    // it->Init(mindist_, &lut_);
    it->ReadCheckpoint(icheck_file);
    n_xlinks_++;
    n_anchors_bound_ += 2;
  }
  /* Close the file */
  icheck_file.close();
}

void CrosslinkManager::InitOutputFiles() {
  if (params_->crosslink.spec_flag) InitSpecFile();
  if (params_->crosslink.checkpoint_flag) InitCheckpoints();
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
