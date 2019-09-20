#include "crosslink.hpp"

Crosslink::Crosslink() : Object() { SetSID(species_id::crosslink); }

void Crosslink::Init(MinimumDistance *mindist, LookupTable *lut) {
  mindist_ = mindist;
  lut_ = lut;
  length_ = -1;
  diameter_ = params_->crosslink.tether_diameter;
  color_ = params_->crosslink.tether_color;
  draw_ = draw_type::_from_string(params_->crosslink.tether_draw_type.c_str());
  rest_length_ = params_->crosslink.rest_length;
  k_on_ = params_->crosslink.k_on;
  k_on_sd_ = params_->crosslink.k_on_sd; // k_on for singly to doubly
  k_off_ = params_->crosslink.k_off;
  k_spring_ = params_->crosslink.k_spring;
  k_align_ = params_->crosslink.k_align;
  f_spring_max_ = params_->crosslink.f_spring_max;
  rcapture_ = params_->crosslink.r_capture;
  fdep_factor_ = params_->crosslink.force_dep_factor;
  polar_affinity_ = params_->crosslink.polar_affinity;
  /* TODO generalize crosslinks to more than two anchors */
  Anchor anchor1;
  Anchor anchor2;
  anchors_.push_back(anchor1);
  anchors_.push_back(anchor2);
  anchors_[0].Init();
  anchors_[1].Init();
  SetSID(species_id::crosslink);
  SetSingly();
}

void Crosslink::UpdatePosition() {
  SetPosition(anchors_[0].GetPosition());
  SetScaledPosition(anchors_[0].GetScaledPosition());
  SetOrientation(anchors_[0].GetOrientation());
  length_ = 0;
}

void Crosslink::GetAnchors(std::vector<Object *> *ixors) {
  if (IsUnbound())
    return;
  std::vector<Object *> ix;
  ix.push_back(&(anchors_[0]));
  if (IsDoubly()) {
    ix.push_back(&(anchors_[1]));
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

void Crosslink::UnbindAnchor(bool second) {
  /* If singly bound, just clear both */
  if (IsSingly()) {
    anchors_[0].Clear();
    anchors_[1].Clear();
    return;
  }
  if (second) {
    anchors_[1].Clear();
  }
  /* When singly bound, only consider first anchor populated */
  else {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
  }
}

/* Get pointer to anchor that is bound. Should only be used for singly bound
 * xlinks. */
Anchor *Crosslink::GetBoundPtr() {
  if (IsDoubly()) {
    warning("GetBoundPtr() called on doubly bound crosslink. This function"
            "should only be called on crosslinks with a free anchor!");
  }
  /* Only first anchor should be bound in the case of singly bound */
  return &(anchors_[0]);
}

/* Perform kinetic monte carlo step of protein with 1 head attached. */
void Crosslink::SinglyKMC() {
  if (!anchors_[0].IsBound()) {
    SetUnbound();
    return;
  }
  /* Must populate filter with 1 for every neighbor, since KMC
  expects a mask. We already guarantee uniqueness, so we won't overcount. */

  double roll = gsl_rng_uniform_pos(rng_.r);
  int head_bound = 0;
  // Set up KMC objects and calculate probabilities
  double unbind_prob = k_off_ * delta_;
  // int n_neighbors = nlist_.size();
  int n_neighbors = neighbors_.NNeighbors();
  std::vector<int> kmc_filter(n_neighbors, 1);
  KMC<Object> kmc_bind(anchors_[0].pos, n_neighbors, rcapture_, delta_, lut_);
  kmc_bind.SetPBCs(n_dim_, space_->n_periodic, space_->unit_cell);
  double kmc_bind_prob = 0;
  std::vector<double> kmc_bind_factor(n_neighbors, k_on_sd_);
  if (n_neighbors > 0) {
    kmc_bind.CalcTotProbsSD(neighbors_.GetNeighborsMem(), kmc_filter,
                            anchors_[0].GetBoundOID(), 0, k_spring_, 1.0,
                            rest_length_, kmc_bind_factor);
    kmc_bind_prob = kmc_bind.getTotProb();
  }
  // Get total probability of changing protein state

  double totProb = kmc_bind_prob + unbind_prob;
  // printf("BindProb: %2.8f\n", kmc_bind_prob);
  // printf("BindProb: %2.2f", kmc_bind_prob;

  int head_activate = -1; // No head activated
  if (totProb > 1.0) {
    // Probability of KMC is greater than one, normalize
    head_activate = (roll < (unbind_prob / totProb)) ? 0 : 1;
    // warning(
    //"Probability of head binding or unbinding in SinglyKMC()"
    //" is >1 (%2.2f). Change time step to prevent this!",
    // totProb);
  } else if (roll < totProb) {
    // Choose which action to perform, bind or unbind
    head_activate = (roll < unbind_prob) ? 0 : 1;
  } else { // No head binds or unbinds
    return;
  }
  // Change status of activated head
  if (head_activate == 0) {
    // Unbind bound head
    anchors_[0].Clear();
    SetUnbound();
  } else if (head_activate == 1) {
    // Bind unbound head
    roll = roll - unbind_prob;
    /* Position on rod where protein will bind with respect to center of rod,
     * passed by reference */
    double bind_lambda;
    /* Pick rod to bind to. Should not give -1 since roll is within the
     * total binding probability. If failure, make sure rolls are shifted
     * properly. */
    int i_bind = kmc_bind.whichRodBindSD(bind_lambda, roll);
    if (i_bind < 0) { // || bind_lambda < 0) {
      printf("i_bind = %d\nbind_lambda = %2.2f\n", i_bind, bind_lambda);
      error_exit("kmc_bind.whichRodBindSD in Crosslink::SinglyKMC"
                 " returned an invalid result!");
    }
    /* KMC returns bind_lambda to be with respect to center of rod. We want it
       to be specified from the tail of the rod to be consistent */
    Object *obj = neighbors_.GetNeighbor(i_bind);
    bind_lambda += 0.5 * obj->GetLength();
    /* KMC can return values that deviate a very small amount from the true rod
       length. Bind to ends if lambda < 0 or lambda > bond_length. */
    if (bind_lambda > obj->GetLength()) {
      bind_lambda = obj->GetLength();
    } else if (bind_lambda < 0) {
      bind_lambda = 0;
    }
    anchors_[1].AttachObjLambda(obj, bind_lambda);
    SetDoubly();
  }
  kmc_filter.clear();
}

/* Perform kinetic monte carlo step of protein with 2 heads of protein
 * object attached. */
void Crosslink::DoublyKMC() {
  if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
    anchors_[0].Clear();
    anchors_[1].Clear();
    SetUnbound();
    return;
  } else if (!anchors_[1].IsBound()) {
    anchors_[1].Clear();
    SetSingly();
    return;
  } else if (!anchors_[0].IsBound()) {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
    SetSingly();
    return;
  }
  /* Check for polar affinities */
  double affinity = 0;
  if (ABS(polar_affinity_) > 1e-8) {
    /* Check whether crosslink is linking two polar or antipolar bonds */
    int po = SIGNOF(dot_product(n_dim_, anchors_[0].GetOrientation(),
                                anchors_[1].GetOrientation()));
    /* polar_affinity_ should range from -1 to 1, with the interval [-1, 0]
     * corresponding to the affinity (between 0 and 1) of binding to antipolar
     * filaments and the interval [0, 1] corresponding to polar aligned
     * filaments. Affinity = 0 means there is no preference, whereas affinity
     * approx 1 means there is a very strong preference for polar filaments and
     * thus a strong likelihood of unbinding from antipolar filaments */
    /* antipolar case */
    if (po < 0 && polar_affinity_ > 0) {
      affinity = polar_affinity_;
    }
    /* polar case */
    else if (po > 0 && polar_affinity_ < 0) {
      affinity = ABS(polar_affinity_);
    }
    if (affinity > 1) {
      affinity = 1;
    }
  }
  /* Calculate force-dependent unbinding for each head */
  double tether_stretch = length_ - rest_length_;
  tether_stretch = (tether_stretch > 0 ? tether_stretch : 0);
  double fdep = fdep_factor_ * 0.5 * k_spring_ * SQR(tether_stretch);
  double totProb;
  if (1 - affinity < 1e-4) {
    totProb = 1;
  } else {
    totProb = 2 * k_off_ * delta_ * exp(fdep) / (1 - affinity);
  }
  double roll = gsl_rng_uniform_pos(rng_.r);
  int head_activate = -1; // No head activated
  if (totProb > 1.0) {    // Probability of KMC is greater than one, normalize
    head_activate = (roll < 0.5) ? 0 : 1;
    // warning(
    //"Probability of head binding or unbinding in DoublyKMC()"
    //" is >1 (%2.2f). Change time step to prevent this!",
    // totProb);
    // printf("tether length: %2.2f\nfdep: %6.2f\n", length_, fdep);
  } else if (roll < totProb) { // Choose which head to unbind
    head_activate = (roll < 0.5) ? 0 : 1;
  } else { // No head unbinds
    return;
  }
  if (head_activate == 0) {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
    SetSingly();
  } else {
    anchors_[1].Clear();
    SetSingly();
  }
}

void Crosslink::CalculateBinding() {
  if (IsSingly()) {
    SinglyKMC();
  } else if (IsDoubly()) {
    DoublyKMC();
  }
  neighbors_.Clear();
}

void Crosslink::ClearNeighbors() { neighbors_.Clear(); }

void Crosslink::UpdateAnchorsToMesh() {
  anchors_[0].UpdateAnchorPositionToMesh();
  anchors_[1].UpdateAnchorPositionToMesh();
}

void Crosslink::UpdateAnchorPositions() {
  anchors_[0].UpdatePosition();
  anchors_[1].UpdatePosition();
}

void Crosslink::ApplyTetherForcesToMesh() {
  anchors_[0].ApplyAnchorForces();
  anchors_[1].ApplyAnchorForces();
}

void Crosslink::UpdateCrosslink() {
  /* Update anchor positions in space to calculate tether forces */
  UpdateAnchorsToMesh();
  /* Check if an anchor became unbound */
  UpdateXlinkState();
  /* If we are doubly-bound, calculate tether forces */
  CalculateTetherForces();
  /* Transfer tether forces to mesh */
  ApplyTetherForcesToMesh();
  /* Have anchors diffuse/walk along mesh */
  UpdateAnchorPositions();
  UpdateXlinkState();
  /* Check for binding/unbinding events using KMC */
  CalculateBinding();
  /* If we are singly bound, update the position of the crosslinker to reflect
   * the singly-bound anchor position for interaction purposes. TODO: Refactor
   * this to use anchors for interactions. */
  if (IsSingly()) {
    UpdatePosition();
  }
}

/* This function ensures that singly-bound crosslinks have anchor[0] bound and
   anchor[1] unbound. */
void Crosslink::UpdateXlinkState() {
  if (IsDoubly() && !anchors_[1].IsBound()) {
    SetSingly();
  } else if (IsDoubly() && !anchors_[0].IsBound()) {
    anchors_[0] = anchors_[1];
    SetSingly();
  }
  if (IsSingly() && !anchors_[0].IsBound()) {
    SetUnbound();
    return;
  }
}

void Crosslink::ZeroForce() {
  std::fill(force_, force_ + 3, 0.0);
  anchors_[0].ZeroForce();
  anchors_[1].ZeroForce();
}

void Crosslink::CalculateTetherForces() {
  ZeroForce();
  if (!IsDoubly())
    return;
  mindist_->ObjectObject(&(anchors_[0]), &(anchors_[1]), &ix);
  /* Check stretch of tether. No penalty for having a stretch < rest_length. ie
   * the spring does not resist compression. */
  length_ = sqrt(ix.dr_mag2);
  double stretch = length_ - rest_length_;
  for (int i = 0; i < params_->n_dim; ++i) {
    orientation_[i] = ix.dr[i] / length_;
    position_[i] = ix.midpoint[i];
  }
  if (stretch > 0) {
    tether_force_ = k_spring_ * stretch;
    for (int i = 0; i < params_->n_dim; ++i) {
      force_[i] = tether_force_ * orientation_[i];
    }
    anchors_[0].AddForce(force_);
    anchors_[1].SubForce(force_);
  }
  UpdatePeriodic();
  /* TODO Apply torques on crosslinks if necessary */
}

/* Attach a crosslink anchor to object in a random fashion. Currently only
 * bonds are used, but can be extended to sites (sphere-like objects) */
void Crosslink::AttachObjRandom(Object *obj) {
  /* Attaching to random bond implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj->GetType() == +obj_type::bond) {
    anchors_[0].AttachObjRandom(obj);
    SetMeshID(obj->GetMeshID());
  } else {
    /* TODO: add binding to sphere or site-like objects */
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
}

void Crosslink::Draw(std::vector<graph_struct *> *graph_array) {
  /* Draw anchors */
  anchors_[0].Draw(graph_array);
  anchors_[1].Draw(graph_array);
  /* Draw tether */
  if (IsDoubly() && length_ > 0) {
    std::copy(scaled_position_, scaled_position_ + 3, g_.r);
    for (int i = space_->n_periodic; i < n_dim_; ++i) {
      g_.r[i] = position_[i];
    }
    // std::copy(position_, position_+3, g_.r);
    std::copy(orientation_, orientation_ + 3, g_.u);
    g_.color = color_;
    if (params_->graph_diameter > 0) {
      g_.diameter = params_->graph_diameter;
    } else {
      g_.diameter = diameter_;
    }
    g_.length = length_;
    g_.draw = draw_;
    graph_array->push_back(&g_);
  }
}

void Crosslink::AddNeighbor(Object *neighbor) {
  neighbors_.AddNeighbor(neighbor);
}

void Crosslink::SetDoubly() {
  state_ = bind_state::doubly;
  SetMeshID(0);
}

void Crosslink::SetSingly() {
  state_ = bind_state::singly;
  SetMeshID(anchors_[0].GetMeshID());
}

void Crosslink::SetUnbound() {
  state_ = bind_state::unbound;
  SetMeshID(0);
}

bool Crosslink::IsDoubly() { return state_ == +bind_state::doubly; }

bool Crosslink::IsSingly() { return state_ == +bind_state::singly; }

bool Crosslink::IsUnbound() { return state_ == +bind_state::unbound; }

void Crosslink::WriteSpec(std::fstream &ospec) {
  if (IsUnbound()) {
    warning("Unbound crosslink tried to WriteSpec!");
    return;
  }
  bool is_doubly = IsDoubly();
  ospec.write(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&length_), sizeof(double));
  double temp[3];
  std::copy(position_, position_ + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  std::copy(orientation_, orientation_ + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  double const *const r0 = anchors_[0].GetPosition();
  std::copy(r0, r0 + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  double const *const u0 = anchors_[0].GetOrientation();
  std::copy(u0, u0 + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  double lambda = anchors_[0].GetBondLambda();
  ospec.write(reinterpret_cast<char *>(&lambda), sizeof(double));
  double const *const r1 = anchors_[1].GetPosition();
  std::copy(r1, r1 + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  double const *const u1 = anchors_[1].GetOrientation();
  std::copy(u1, u1 + 3, temp);
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  lambda = anchors_[1].GetBondLambda();
  ospec.write(reinterpret_cast<char *>(&lambda), sizeof(double));
}

void Crosslink::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  bool is_doubly;
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length_), sizeof(double));
  double temp[3];
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  SetPosition(temp);
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  SetOrientation(temp);
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  anchors_[0].SetPosition(temp);
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  anchors_[0].SetOrientation(temp);
  double lambda;
  ispec.read(reinterpret_cast<char *>(&lambda), sizeof(double));
  anchors_[0].SetBondLambda(lambda);
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  anchors_[1].SetPosition(temp);
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&temp[i]), sizeof(double));
  }
  anchors_[1].SetOrientation(temp);
  ispec.read(reinterpret_cast<char *>(&lambda), sizeof(double));
  anchors_[1].SetBondLambda(lambda);
  anchors_[0].UpdatePeriodic();
  anchors_[1].UpdatePeriodic();
  anchors_[0].SetBound();
  if (is_doubly) {
    SetDoubly();
    UpdatePeriodic();
    anchors_[1].SetBound();
  } else {
    UpdatePosition();
  }
}

void Crosslink::WriteCheckpoint(std::fstream &ocheck) {
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char *>(rng_state), rng_size);
  WriteSpec(ocheck);
}

void Crosslink::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof())
    return;
  void *rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char *>(rng_state), rng_size);
  ReadSpec(icheck);
}
