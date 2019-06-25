#include "crosslink.h"

Crosslink::Crosslink() : Object() {
  SetSID(species_id::crosslink);
}

void Crosslink::Init(MinimumDistance * mindist, LookupTable * lut) {
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
  /* TODO generalize crosslinks to more than two anchors */
  Anchor anchor1, anchor2;
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
Anchor * Crosslink::GetBoundPtr() {
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
    anchors_[0].Clear();
    SetUnbound();
    return;
  }
  double roll = gsl_rng_uniform_pos(rng_.r);
  int head_bound = 0;
  // Set up KMC objects and calculate probabilities
  double unbind_prob = k_off_ * delta_;
  int n_neighbors = nlist_.size();
  //for (int i=0; i<n_dim_; ++i) {
    //printf("%2.2f, ", anchors_[0].GetPosition()[i]);
  //}
  //printf("\n");
  //for (int i=0; i<n_dim_; ++i) {
    //printf("%2.2f, ", anchors_[0].pos[i]);
  //}
  //printf("\n");
  //printf("%d\n", anchors_[0].GetBoundOID());
  KMC<Object> kmc_bind(anchors_[0].pos, nlist_.size(), 
                     rcapture_, lut_, n_dim_,
                     space_->n_periodic, space_->unit_cell);
  //KMC<Object> kmc_bind(anchors_[0].GetPosition(), nlist_.size(), rcapture_);
  double kmc_bind_prob = 0;
  if (n_neighbors > 0) {
    kmc_bind.CalcTotProbsSD(&(nlist_[0]), kmc_filter_, anchors_[0].GetBoundOID(),
                            0, k_spring_, 1.0, rest_length_,
                            k_on_sd_*delta_);
    kmc_bind_prob = kmc_bind.getTotProb();
  }
  // Get total probability of changing protein statek

  double totProb = kmc_bind_prob + unbind_prob;
  //printf("BindProb: %2.8f\n", kmc_bind_prob);
  //printf("BindProb: %2.2f", kmc_bind_prob;

  int head_activate = -1; // No head activated
  if (totProb > 1.0) {    
    // Probability of KMC is greater than one, normalize
    head_activate =  (roll < (unbind_prob / totProb)) ? 0 : 1;
    warning("Probability of head binding or unbinding in SinglyKMC()"
      " is >1 (%2.2f). Change time step to prevent this!", totProb);
  } 
  else if (roll < totProb) { 
    // Choose which action to perform, bind or unbind
    head_activate = (roll < unbind_prob) ? 0 : 1;
  } 
  else { // No head binds or unbinds
    return;
  }
  // Change status of activated head
  if (head_activate == 0) {
    // Unbind bound head
    anchors_[0].Clear();
    SetUnbound();
  } 
  else if (head_activate == 1) { 
    // Bind unbound head
    roll = roll - unbind_prob;
    /* Position on rod where protein will bind, passed by reference */
    double bind_lambda;
    /* Pick rod to bind to. Should not give -1 since roll is within the 
     * total binding probability. If failure, make sure rolls are shifted 
     * properly. */
    int i_bind = kmc_bind.whichRodBindSD(bind_lambda, roll);
    if (i_bind < 0) {// || bind_lambda < 0) {
      printf("i_bind = %d\nbind_lambda = %2.2f\n", i_bind, bind_lambda);
      error_exit("kmc_bind.whichRodBindSD in Crosslink::SinglyKMC"
          " returned an invalid result!");
    }
    anchors_[1].AttachObjLambda(nlist_[i_bind], bind_lambda);
    SetDoubly();
  }
}

/* Perform kinetic monte carlo step of protein with 2 heads of protein 
 * object attached. */
void Crosslink::DoublyKMC() {
  if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
    anchors_[0].Clear();
    anchors_[1].Clear();
    SetUnbound();
    return;
  }
  else if (!anchors_[1].IsBound()) {
    anchors_[1].Clear();
    SetSingly();
    return;
  }
  else if (!anchors_[0].IsBound()) {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
    SetSingly();
    return;
  }
  /* Calculate force-dependent unbinding for each head */
  double tether_stretch = length_ - rest_length_;
  tether_stretch = (tether_stretch > 0 ? tether_stretch : 0);
  double fdep = fdep_factor_ * 0.5 * k_spring_ * SQR(tether_stretch);
  double totProb = 2 * k_off_ * delta_ * exp(fdep);
  double roll = gsl_rng_uniform_pos(rng_.r);
  int head_activate = -1; // No head activated
  if (totProb > 1.0) {    // Probability of KMC is greater than one, normalize
    head_activate = (roll < 0.5) ? 0 : 1;
    warning("Probability of head binding or unbinding in DoublyKMC()"
      " is >1 (%2.2f). Change time step to prevent this!", totProb);
    printf("tether length: %2.2f\nfdep: %6.2f\n", length_, fdep);
  } 
  else if (roll < totProb) { // Choose which head to unbind
    head_activate = (roll < 0.5) ? 0 : 1;
  } 
  else { // No head unbinds
    return;
  }
  if (head_activate == 0) {
    anchors_[0] = anchors_[1];
    anchors_[1].Clear();
    SetSingly();
  }
  else {
    anchors_[1].Clear();
    SetSingly();
  }
}

void Crosslink::CalculateBinding() {
  if (IsSingly()) {
    SinglyKMC();
  }
  else if (IsDoubly()) {
    DoublyKMC();
  }
  nlist_.clear();
  kmc_filter_.clear();
}
void Crosslink::CalcBinding() {
  if (IsDoubly()) {
    //[> If our second anchor became unbound <]
    if (length_ > rcapture_ || tether_force_ > f_spring_max_) {
      //[> Designate the survivor of the crosslink breaking <]
      int which = gsl_rng_uniform_int(rng_.r, 2);
      anchors_[0] = anchors_[which];
      anchors_[1].Clear();
      SetSingly();
      return;
    }
    if (!anchors_[1].IsBound()) {
      anchors_[1].Clear();
      SetSingly();
    }
    //[> If our first anchor became unbound <]
    if (!anchors_[0].IsBound()) {
      //[> Swap anchors so that anchor[0] remains the one singly bound
       //* anchor. If anchor 2 became unbound also, we'll
         //catch that in a moment */
      anchors_[0] = anchors_[1];
      anchors_[1].Clear();
      SetSingly();
    }
    //[> If both anchors became unbound <]
    if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
      SetUnbound();
    }
  }
  //[> Check if our singly-bound anchor became unbound <]
  else if (!anchors_[0].IsBound()) {
    SetUnbound();
    anchors_[0].Clear();
  }
  else if (IsSingly()) {
    AttemptCrosslink();
  }
  nlist_.clear();
  kmc_filter_.clear();
}

void Crosslink::UpdateCrosslink() {
  /* If we are doubly-bound, calculate tether forces */
  if (IsDoubly()) {
    CalculateTetherForces();
  }
  anchors_[0].UpdatePosition();
  anchors_[1].UpdatePosition();
  /* Check if any of our doubly-bound anchors became unbound */
  CalculateBinding(); // Using KMC
  //CalcBinding(); // Using naive binding rules
  /* If we are singly bound, update the position of the crosslinker
   * to reflect the singly-bound anchor position for interaction
   * purposes */
  if (IsSingly()) {
    UpdatePosition();
  }
}

void Crosslink::AttemptCrosslink() {
  int n_neighbors = nlist_.size();
  if (n_neighbors == 0) {
    return;
  }
  /* Do stupidest possible thing, and bind to whatever with some probability
   * and in a random fashion */
  double p_bind = k_on_*delta_;
  if (gsl_rng_uniform_pos(rng_.r) < p_bind) {
    int i_bind = gsl_rng_uniform_int(rng_.r, n_neighbors);
    Object * obj = nlist_[i_bind];
    mindist_->ObjectObject(&(anchors_[0]), obj, &ix);
    if (ix.dr_mag2 > SQR(rcapture_)) {
      return;
    }
    if (obj->GetType() == +obj_type::bond) {
      anchors_[1].AttachObjRandom(obj);
    }
    else {
      /* TODO: add binding to sphere or site-like objects */
      error_exit("Crosslink binding to non-bond objects not yet implemented.");
    }
    SetDoubly();
  }
}

void Crosslink::CalculateTetherForces() {
  mindist_->ObjectObject(&(anchors_[0]), &(anchors_[1]), 
      &ix);
  /* Check stretch of tether. No penalty for having a stretch < rest_length. ie
   * the spring does not resist compression. */
  length_ = sqrt(ix.dr_mag2);
  double stretch = length_ - rest_length_;
  std::fill(force_, force_+3, 0.0);
  for (int i=0; i<params_->n_dim; ++i) {
    orientation_[i] = ix.dr[i]/length_;
    position_[i] = ix.midpoint[i];
  }
  if (stretch > 0) {
    tether_force_ = k_spring_ * stretch;
    for (int i=0; i<params_->n_dim; ++i) {
      force_[i] = tether_force_ * orientation_[i];
    }
    anchors_[0].AddForce(force_);
    anchors_[1].SubForce(force_);
    anchors_[0].ApplyAnchorForces();
    anchors_[1].ApplyAnchorForces();
  }
  UpdatePeriodic();
  /* TODO Apply torques on crosslinks if necessary */
}

/* Attach a crosslink anchor to object in a random fashion. Currently only
 * bonds are used, but can be extended to sites (sphere-like objects) */
void Crosslink::AttachObjRandom(Object * obj) {
  /* Attaching to random bond implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj->GetType() == +obj_type::bond) {
    anchors_[0].AttachObjRandom(obj);
    SetMeshID(obj->GetMeshID());
  }
  else {
    /* TODO: add binding to sphere or site-like objects */
    error_exit("Crosslink binding to non-bond objects not yet implemented.");
  }
}

void Crosslink::Draw(std::vector<graph_struct*> * graph_array) {
  /* Draw anchors */
  anchors_[0].Draw(graph_array);
  anchors_[1].Draw(graph_array);
  /* Draw tether */
  if (IsDoubly() && length_ > 0) {
    std::copy(scaled_position_,scaled_position_+3, g_.r);
    for (int i=space_->n_periodic; i<n_dim_; ++i) {
      g_.r[i] = position_[i];
    }
    //std::copy(position_, position_+3, g_.r);
    std::copy(orientation_, orientation_+3, g_.u);
    g_.color = color_;
    if (params_->graph_diameter > 0) {
      g_.diameter = params_->graph_diameter;
    }
    else {
      g_.diameter = diameter_;
    }
    g_.length = length_;
    g_.draw = draw_;
    graph_array->push_back(&g_);
  }
}

void Crosslink::AddNeighbor(Object * neighbor) {
  /* TODO Prevent racy conditions */
  //std::lock_guard<std::mutex> lk(xlink_mtx_);
  nlist_.push_back(neighbor);
  /* Must populate filter with 1 for every neighbor, since KMC
  expects a mask. We already guarantee uniqueness, so we won't overcount. */
  kmc_filter_.push_back(1);
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

bool Crosslink::IsDoubly() {
  return state_ == +bind_state::doubly;
}

bool Crosslink::IsSingly() {
  return state_ == +bind_state::singly;
}

bool Crosslink::IsUnbound() {
  return state_ == +bind_state::unbound;
}

void Crosslink::WriteSpec(std::fstream &ospec) {
  if (IsUnbound()) {
    warning("Unbound crosslink tried to WriteSpec!");
    return;
  }
  bool is_doubly = IsDoubly();
  ospec.write(reinterpret_cast<char*>(&is_doubly), sizeof(bool));
  ospec.write(reinterpret_cast<char*>(&diameter_), sizeof(double));
  ospec.write(reinterpret_cast<char*>(&length_), sizeof(double));
  double temp[3];
  std::copy(position_, position_+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  std::copy(orientation_, orientation_+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  double const * const r0 = anchors_[0].GetPosition();
  std::copy(r0, r0+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  double const * const u0 = anchors_[0].GetOrientation();
  std::copy(u0, u0+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  double lambda = anchors_[0].GetBondLambda();
  ospec.write(reinterpret_cast<char*>(&lambda), sizeof(double));
  double const * const r1 = anchors_[1].GetPosition();
  std::copy(r1, r1+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  double const * const u1 = anchors_[1].GetOrientation();
  std::copy(u1, u1+3, temp);
  for (int i=0; i<3; ++i) {
    ospec.write(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  lambda = anchors_[1].GetBondLambda();
  ospec.write(reinterpret_cast<char*>(&lambda), sizeof(double));
}

void Crosslink::ReadSpec(std::fstream &ispec) {
  if (ispec.eof()) return;
  bool is_doubly;
  ispec.read(reinterpret_cast<char*>(&is_doubly), sizeof(bool));
  ispec.read(reinterpret_cast<char*>(&diameter_), sizeof(double));
  ispec.read(reinterpret_cast<char*>(&length_), sizeof(double));
  double temp[3];
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  SetPosition(temp);
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  SetOrientation(temp);
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  anchors_[0].SetPosition(temp);
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  anchors_[0].SetOrientation(temp);
  double lambda;
  ispec.read(reinterpret_cast<char*>(&lambda), sizeof(double));
  anchors_[0].SetBondLambda(lambda);
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  anchors_[1].SetPosition(temp);
  for (int i=0; i<3; ++i) {
    ispec.read(reinterpret_cast<char*>(&temp[i]), sizeof(double));
  }
  anchors_[1].SetOrientation(temp);
  ispec.read(reinterpret_cast<char*>(&lambda), sizeof(double));
  anchors_[1].SetBondLambda(lambda);
  anchors_[0].UpdatePeriodic();
  anchors_[1].UpdatePeriodic();
  anchors_[0].SetBound();
  if (is_doubly) {
    SetDoubly();
    UpdatePeriodic();
    anchors_[1].SetBound();
  }
  else {
    UpdatePosition();
  }
}

void Crosslink::WriteCheckpoint(std::fstream &ocheck) {
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char*>(rng_state), rng_size);
  WriteSpec(ocheck);
}

void Crosslink::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char*>(rng_state), rng_size);
  ReadSpec(icheck);
}

