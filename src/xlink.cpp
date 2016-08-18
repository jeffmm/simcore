#include "xlink.h"

#include <iomanip>

void Xlink::Init() {
  Composite::InsertRandom(length_+diameter_);
  // Give each bead the location of the main position!
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->SetPosition(position_);
    i_bead->SetPrevPosition(position_);
    i_bead->SetDiameter(diameter_);
    i_bead->SetDiffusion();
    i_bead->SetSpace(space_);
    i_bead->UpdatePeriodic();
    i_bead->SetBound(false);
  }

  UpdateOrientation();
}

void Xlink::InitConfigurator(const double* const x, const double diameter) {
  SetPosition(x);
  diameter_ = diameter;
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->SetPosition(position_);
    i_bead->SetPrevPosition(position_);
    i_bead->SetDiameter(diameter_);
    i_bead->SetDiffusion();
    i_bead->SetSpace(space_);
    i_bead->UpdatePeriodic();
    i_bead->SetBound(false);
  }

  UpdateOrientation();
}

void Xlink::UpdateOrientation() {
  double const * const r1=elements_[0].GetPosition();
  double const * const r2=elements_[1].GetPosition();
  double const * const s1=elements_[0].GetScaledPosition();
  double const * const s2=elements_[1].GetScaledPosition();
  length_ = 0;
  double dr[3];
  separation_vector(n_dim_, space_->n_periodic, r1, s1, r2, s2, space_->unit_cell, dr);
  for (int i=0; i<n_dim_; ++i) {
    length_ += SQR(dr[i]);
  }
  if (length_ <= 0.0)
    length_ = 1.0;
  length_=sqrt(length_);
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = (dr[i])/length_;
    position_[i] = r1[i] + 0.5*length_*orientation_[i];
  }
  UpdatePeriodic();
  double orientation_loc[3];
  orientation_loc[0] = 0;
  orientation_loc[1] = 1;
  orientation_loc[2] = 1;
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->SetOrientation(orientation_loc); 
  }
}

void Xlink::CheckBoundState() {
  bound_ = unbound;
  auto head0 = elements_.begin();
  auto head1 = elements_.begin()+1;
  bool bound0 = head0->GetBound();
  bool bound1 = head1->GetBound();
  if (!bound0 && !bound1) {
    uinternal_ = 0.0;
  }
  if (bound0 && !bound1) {
    bound_ = singly;
    uinternal_ = 0.0;
  }
  if (!bound0 && bound1) {
    bound_ = singly;
    uinternal_ = 0.0;
  }
  if (bound0 && bound1) {
    bound_ = doubly;
  }
}

std::pair<bool,bool> Xlink::GetBoundHeads(XlinkHead **freehead, XlinkHead **boundhead) {
  std::pair<bool,bool> nbound;
  auto head0 = elements_.begin();
  auto head1 = elements_.begin()+1;
  if (!head0->GetBound() && !head1->GetBound()) {
    // Neither is bound, just return
    nbound.first = false;
    nbound.second = false;
  } else if (head0->GetBound() && !head1->GetBound()) {
    // Head 0 is bound and not head1
    nbound.first = true;
    nbound.second = false;
    (*freehead) = &(*head1);
    (*boundhead) = &(*head0);
  } else if (!head0->GetBound() && head1->GetBound()) {
    nbound.first = false;
    nbound.second = true;
    (*freehead) = &(*head0);
    (*boundhead) = &(*head1);
  } else {
    // Both heads are bound
    nbound.first = true;
    nbound.second = true;
  }
  return nbound;
}

void Xlink::UpdatePositionMP() {
  CheckBoundState();
  Integrate();
}

void Xlink::ApplyInteractions() {
   //We have no internal constraints
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->ApplyInteractions();
  }
}

void Xlink::Integrate() {
  // Different actions due to being free, singly, or doubly bound
  switch(bound_) {
    case unbound:
      DiffuseXlink();
      break;
  }
  /*double pos[3] = {0, 0, 0};
  // the beads know how to update position and periodicity
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead) {
    i_bead->UpdatePositionMP();
  }*/
  UpdateOrientation();
}

void Xlink::DiffuseXlink() {
  // Diffuse based on the first head
  auto head0 = elements_.begin();
  auto head1 = elements_.begin()+1;
  double oldpos[3];
  std::copy(head0->GetRigidPosition(), head0->GetRigidPosition()+n_dim_, oldpos);
  head0->UpdatePositionMP();
  head0->SetPrevPosition(oldpos);
  head1->SetPosition(head0->GetRigidPosition());
  head1->SetPrevPosition(oldpos);
  head1->AddDr();
  SetPosition(head0->GetRigidPosition());
  SetPrevPosition(oldpos);
  head0->UpdatePeriodic();
  head1->UpdatePeriodic();
  UpdatePeriodic();
}

void Xlink::Draw(std::vector<graph_struct*> * graph_array) {
  // draw half of rod coming from each bead (looks good for periodic boundaries)
  double const * const r1 = elements_[0].GetPosition();
  double const * const r2 = elements_[1].GetPosition();
  double r[3];
  for (int i=0; i<n_dim_; ++i)
    r[i] = r1[i] + 0.25*length_*orientation_[i];
  std::copy(r, r+3, g_.r);
  std::copy(orientation_, orientation_+3, g_.u);
  g_.length = 0.5 * length_;
  g_.diameter = 0.2*diameter_;
  graph_array->push_back(&g_);
  for (int i=0; i<n_dim_; ++i)
    r[i] = r2[i] - 0.25*length_*orientation_[i];
  std::copy(r, r+3, g2_.r);
  std::copy(orientation_, orientation_+3, g2_.u);
  g2_.length = 0.5 * length_;
  g2_.diameter = 0.2*diameter_;
  graph_array->push_back(&g2_);
  for (auto i_bead = elements_.begin(); i_bead != elements_.end(); ++i_bead)
    i_bead->Draw(graph_array);
}


// Species specifics
void XlinkSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "Xlink species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";
  std::string insertion_type;
  insertion_type = node["xlink"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["xlink"]["properties"]["overlap"].as<bool>();
  std::cout << "   can overlap:    " << (can_overlap ? "true" : "false") << std::endl;

  if (insertion_type.compare("xyz") == 0) {
    if (!can_overlap) {
      std::cout << "Warning, location insertion overrides overlap\n";
      can_overlap = true;
    }
    int nxlinks = (int)node["xlink"]["xit"].size();
    std::cout << "   nxlinks: " << nxlinks << std::endl;
    params_->n_xlink = nxlinks;
    for (int ix = 0; ix < nxlinks; ++ix) {
      double x[3] = {0.0, 0.0, 0.0};
      double diameter = 0.0;
      x[0] = node["xlink"]["xit"][ix]["x"][0].as<double>();
      x[1] = node["xlink"]["xit"][ix]["x"][1].as<double>();
      x[2] = node["xlink"]["xit"][ix]["x"][2].as<double>();
      std::cout << "   x(" << x[0] << ", " << x[1] << ", " << x[2] << ")\n";
      diameter = node["xlink"]["xit"][ix]["diameter"].as<double>();
      std::cout << "   diameter[" << diameter << "]\n";

      Xlink *member = new Xlink(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(x, diameter);
      member->Dump();
      members_.push_back(member);
    }
  } else if (insertion_type.compare("random") == 0) {
    int nxlinks     = node["xlink"]["xit"]["num"].as<int>();
    double diameter = node["xlink"]["xit"]["diameter"].as<double>();

    std::cout << std::setw(25) << std::left << "   n xlinks:" << std::setw(10)
      << std::left << nxlinks << std::endl;
    std::cout << std::setw(25) << std::left << "   diameter:" << std::setw(10)
      << std::left << diameter << std::endl;

    params_->n_xlink = nxlinks;
    params_->xlink_diameter = diameter;

    if (!can_overlap) {
      error_exit("ERROR, xlinks always allowed to overlap\n");
    }
    for (int i = 0; i < nxlinks; ++i) {
      Xlink* member = new Xlink(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      members_.push_back(member);
    }
  } else {
    printf("nope, not yet\n");
    exit(1);
  }
}
