#include "bond.hpp"

/**************************
** Bond member functions **
**************************/
Bond::Bond(unsigned long seed) : Object(seed) { type_ = obj_type::bond; }

void Bond::Init(Site *s1, Site *s2) {
  s1->AddBond(this, OUTGOING);
  s2->AddBond(this, INCOMING);
  sites_[0] = s1;
  sites_[1] = s2;
  ReInit();
  SetEquilLength(length_);
  Logger::Trace("Initializing bond at [%2.2f %2.2f %2.2f] with orientation "
                "[%2.2f %2.2f %2.2f] and length %2.2f",
                position_[0], position_[1], position_[2], orientation_[0],
                orientation_[1], orientation_[2], length_);
}

void Bond::ReInit() {
  double const *const r1 = sites_[0]->GetPosition();
  double const *const r2 = sites_[1]->GetPosition();
  diameter_ = sites_[0]->GetDiameter();
  length_ = 0;
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] = r2[i] - r1[i];
    length_ += orientation_[i] * orientation_[i];
  }
  length_ = sqrt(length_);
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = r1[i] + 0.5 * orientation_[i];
    orientation_[i] /= length_;
  }
  UpdatePeriodic();
}

void Bond::SetEquilLength(double el) {
  equil_length_ = el;
}

void Bond::SetMeshPtr(Object *obj_ptr) { mesh_ptr_ = obj_ptr; }

void Bond::SetBondNumber(int bnum) { bond_number_ = bnum; }

int const Bond::GetBondNumber() { return bond_number_; }

Site *Bond::GetSite(int i) {
  if (i < 0 || i > 1) {
    Logger::Error("Requested adjacent site out of bounds!");
  }
  return sites_[i];
}
Bond *Bond::GetNeighborBond(int i) {
  if (i < 0 || i > 1) {
    Logger::Error("Requested neighboring bond out of bounds!");
  }
  return sites_[i]->GetOtherBond(GetOID());
}

directed_bond Bond::GetNeighborDirectedBond(int i) {
  if (i < 0 || i > 1) {
    Logger::Error("Requested neighboring bond out of bounds!");
  }
  return sites_[i]->GetOtherDirectedBond(GetOID());
}

void Bond::Report() {
  fprintf(stderr, "  Bond:\n");
  Object::Report();
}

void Bond::ReportSites() {
  Report();
  fprintf(stderr, "    Reporting sites:\n");
  for (int i = 0; i < 2; ++i) {
    sites_[i]->Report();
  }
}

void Bond::Draw(std::vector<graph_struct *> &graph_array) {
  std::copy(scaled_position_, scaled_position_ + 3, g_.r);
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_ + 3, g_.u);
  g_.color = color_;
  if (params_->graph_diameter > 0) {
    g_.diameter = params_->graph_diameter;
  } else {
    g_.diameter = diameter_;
  }
  g_.length = length_;
  g_.draw = draw_;
  if (has_overlap_ && params_->highlight_overlaps) {
    g_.draw = draw_type::bw;
    g_.diameter = 2 * diameter_;
  }
  int flock_type = GetFlockType();
  if (flock_type && params_->highlight_flock) {
    g_.draw = draw_type::fixed;
    if (flock_type == 1) {
      // Part of flock interior
      g_.color = params_->flock_color_int;
    } else if (flock_type == 2) {
      // Part of flock exterior
      g_.color = params_->flock_color_ext;
    } else {
      Logger::Warning("Unexpected flock parameter value in Bond::Draw");
    }
    g_.diameter = 2 * diameter_;
    SetFlockType(0);
  }
  graph_array.push_back(&g_);
  HasOverlap(false);
}

bool Bond::HasNeighbor(int other_oid) {
  return (sites_[0]->HasNeighbor(other_oid) ||
          sites_[1]->HasNeighbor(other_oid));
}

double const Bond::GetOrientationCorrelation() {
  return dot_product(n_dim_, orientation_, orientation_0_);
}

void Bond::ZeroOrientationCorrelation() {
  std::copy(orientation_, orientation_ + 3, orientation_0_);
}

