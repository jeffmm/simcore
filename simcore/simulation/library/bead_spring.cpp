#include "bead_spring.hpp"

BeadSpring::BeadSpring() : Mesh() {}

void BeadSpring::SetParameters() {
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  persistence_length_ = sparams_->persistence_length;
  diameter_ = sparams_->diameter;
  max_bond_length_ = sparams_->max_bond_length;
  bond_rest_length_ = sparams_->bond_rest_length;
  bond_spring_ = sparams_->bond_spring;
  driving_factor_ = sparams_->driving_factor;
  stoch_flag_ =
      params_->stoch_flag;  // determines whether we are using stochastic forces
  eq_steps_ = params_->n_steps_equil;
  eq_steps_count_ = 0;
}

void BeadSpring::Init(bead_spring_parameters sparams) {
  sparams_ = sparams;
  SetParameters();
  InitElements();
  InsertBeadSpring();
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
}

void BeadSpring::InitElements() {
  n_bonds_max_ = (int)ceil(length_ / (bond_rest_length_));
  Reserve(n_bonds_max_);
  // Allocate control structures
  int n_sites_max = n_bonds_max_ + 1;
  cos_thetas_.resize(n_sites_max - 2);  // max_sites-2
}

void BeadSpring::GenerateProbableOrientation() {
  /* This updates the current orientation with a generated probable orientation
  where we generate random theta pulled from probability distribution P(th) =
  exp(k cos(th)) where k is the persistence_length of the filament If k is too
  large, there is enormous imprecision in this calculation since sinh(k) is very
  large so to fix this I introduce an approximate distribution that is valid
  for large k */
  double theta;
  if (persistence_length_ == 0) {
    theta = rng_.RandomUniform() * M_PI;
  } else if (persistence_length_ < 100) {
    theta = acos(log(exp(-persistence_length_ / bond_length_) +
                     2.0 * rng_.RandomUniform() *
                         sinh(persistence_length_ / bond_length_)) /
                 (persistence_length_ / bond_length_));
  } else {
    theta = acos((log(2.0 * rng_.RandomUniform()) - log(2.0) +
                  persistence_length_ / bond_length_) /
                 (persistence_length_ / bond_length_));
  }
  double new_orientation[3] = {0, 0, 0};
  if (n_dim_ == 2) {
    theta = (rng_.RandomInt(2) == 0 ? -1 : 1) * theta;
    new_orientation[0] = cos(theta);
    new_orientation[1] = sin(theta);
  } else {
    double phi = rng_.RandomUniform() * 2.0 * M_PI;
    new_orientation[0] = sin(theta) * cos(phi);
    new_orientation[1] = sin(theta) * sin(phi);
    new_orientation[2] = cos(theta);
  }
  rotate_vector_relative(n_dim_, new_orientation, orientation_);
  std::copy(new_orientation, new_orientation + 3, orientation_);
}

void BeadSpring::InsertFirstBond() {
  if (sparams_->insertion_type.compare("random") == 0) {
    InitRandomSite(diameter_);
    AddRandomBondToTip(bond_length_);
  } else if (sparams_->insertion_type.compare("random_oriented") ==
             0) {
    InitRandomSite(diameter_);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = (rng_.RandomUniform() > 0.5 ? 1.0 : -1.0);
    AddBondToTip(orientation_, bond_length_);
  } else if (sparams_->insertion_type.compare("centered_oriented") ==
             0) {
    std::fill(position_, position_ + 3, 0.0);
    position_[n_dim_ - 1] = -0.5 * length_;
    InitSiteAt(position_, diameter_);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
    AddBondToTip(orientation_, bond_length_);
  } else if (sparams_->insertion_type.compare("centered_random") ==
             0) {
    generate_random_unit_vector(n_dim_, orientation_, rng_.r);
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] = -0.5 * length_ * orientation_[i];
    }
    InitSiteAt(position_, diameter_);
    AddBondToTip(orientation_, bond_length_);
  }
}

void BeadSpring::InsertBeadSpring() {
  bond_length_ = bond_rest_length_;
  Clear();
  InsertFirstBond();
  SetOrientation(bonds_[n_bonds_ - 1].GetOrientation());
  bool probable_orientation =
      (sparams_->insertion_type.compare("simple_crystal") != 0 &&
       sparams_->insertion_type.compare("random_oriented") != 0);
  for (int i = 0; i < n_bonds_max_ - 1; ++i) {
    if (probable_orientation) {
      GenerateProbableOrientation();
    }
    AddBondToTip(orientation_, bond_length_);
  }
  for (bond_iterator bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    bond->SetColor(color_, draw_);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
}

void BeadSpring::InsertAt(double *new_pos, double *u) {
  bond_length_ = bond_rest_length_;
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = new_pos[i] - 0.5 * length_ * u[i];
  }
  InitSiteAt(position_, diameter_);
  std::copy(u, u + 3, orientation_);
  AddBondToTip(orientation_, bond_length_);
  SetOrientation(bonds_[n_bonds_ - 1].GetOrientation());
  for (int i = 0; i < n_bonds_max_ - 1; ++i) {
    if (sparams_->insertion_type.compare("simple_crystal") != 0) {
      GenerateProbableOrientation();
    }
    AddBondToTip(orientation_, bond_length_);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
}

void BeadSpring::SetDiffusion() {
  rand_sigma_ = sqrt(24.0 * diameter_ / delta_);
}

double const BeadSpring::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_ * length_ + 0.25 * M_PI * diameter_ * diameter_;
  }
  if (n_dim_ == 3) {
    return 0.25 * M_PI * diameter_ * diameter_ * length_ +
           1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_;
  }
}

void BeadSpring::UpdatePosition(bool midstep) {
  // ApplyForcesTorques();
  Integrate(midstep);
  eq_steps_count_++;
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
void BeadSpring::Integrate(bool midstep) {
  CalculateAngles();
  CalculateTangents();
  if (midstep) {
    CalcRandomForces();
    UpdatePrevPositions();
  }
  AddRandomForces();
  AddBendingForces();
  AddSpringForces();
  UpdateSitePositions(midstep);
  UpdateBondPositions();
  UpdateSiteOrientations();
}

// void BeadSpring::UpdatePrevPositions() {
// for (int i_site=0;i_site<n_sites_;++i_site) {
// sites_[i_site].SetPrevPosition(sites_[i_site].GetPosition());
//}
//}

void BeadSpring::UpdateSiteOrientations() {
  for (int i = 0; i < n_sites_ - 1; ++i) {
    sites_[i].SetOrientation(bonds_[i].GetOrientation());
  }
  sites_[n_sites_ - 1].SetOrientation(bonds_[n_bonds_ - 1].GetOrientation());
}

void BeadSpring::CalculateAngles() {
  double cos_angle;
  for (int i_site = 0; i_site < n_sites_ - 2; ++i_site) {
    double const *const u1 = sites_[i_site].GetOrientation();
    double const *const u2 = sites_[i_site + 1].GetOrientation();
    cos_angle = dot_product(n_dim_, u1, u2);
    cos_thetas_[i_site] = cos_angle;
  }
}

void BeadSpring::CalculateTangents() {
  for (auto it = sites_.begin(); it != sites_.end(); ++it) {
    it->CalcTangent();
  }
}

void BeadSpring::CalcRandomForces() {
  if (!stoch_flag_) return;
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    for (int i = 0; i < n_dim_; ++i) {
      double kick = rng_.RandomUniform() - 0.5;
      force_[i] = kick * rand_sigma_;
    }
    sites_[i_site].SetRandomForce(force_);
  }
}

void BeadSpring::AddRandomForces() {
  if (!stoch_flag_) return;
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    sites_[i_site].AddRandomForce();
  }
}

void BeadSpring::AddBendingForces() {
  /* These forces are calculated from the potential U = k(theta-pi)^2 where
   * theta is the angle between two neighboring bonds consisting of three
   * beads.  I will include the derivation of the site forces in the
   * documentation (JMM 20171215) */
  for (int i_theta = 0; i_theta < n_bonds_ - 1; ++i_theta) {
    double udotu = cos_thetas_[i_theta];
    double kappa;
    if (1 - udotu < 1e-10) {
      kappa = 0;
    } else {
      kappa = -2 * persistence_length_ * (acos(cos_thetas_[i_theta]) - M_PI) /
              sqrt(1 - SQR(cos_thetas_[i_theta]));
    }
    double linv1 = 1.0 / bonds_[i_theta].GetLength();
    double linv2 = 1.0 / bonds_[i_theta + 1].GetLength();
    double const *const u1 = sites_[i_theta].GetOrientation();
    double const *const u2 = sites_[i_theta + 1].GetOrientation();
    // Forces for the first of three sites, which is part of the force on the
    // second site
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] = kappa * linv1 * linv1 * (-u2[i] + udotu * u1[i]);
    }
    sites_[i_theta].AddForce(force_);
    sites_[i_theta + 1].SubForce(force_);
    // Forces for the third of three sites, which is part of the force on the
    // second site
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] = kappa * linv2 * linv2 * (u1[i] - udotu * u2[i]);
    }
    sites_[i_theta + 2].AddForce(force_);
    sites_[i_theta + 1].SubForce(force_);
  }
}

void BeadSpring::AddSpringForces() {
  /* FENE springs have the potential U(r) = -1/2 k r_max^2 ln(1 -
   * ((r-r0)/r_max)^2) */
  for (int i_bond = 0; i_bond < n_bonds_; ++i_bond) {
    double r_diff = bonds_[i_bond].GetLength() - bond_rest_length_;
    double f_mag =
        bond_spring_ * r_diff / (1 - SQR(r_diff) / SQR(max_bond_length_));
    // double f_mag = bond_spring_ * (r-max_bond_length_);
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] = f_mag * sites_[i_bond].GetOrientation()[i];
    }
    sites_[i_bond + 1].SubForce(force_);
    sites_[i_bond].AddForce(force_);
  }

  // for(int i_bond=0; i_bond<n_bonds_; ++i_bond) {
  // double r = bonds_[i_bond].GetLength();
  // double f_mag = bond_spring_ * r / (1 - SQR(r)/SQR(max_bond_length_));
  // for (int i=0;i<n_dim_;++i) {
  // force_[i] = f_mag * sites_[i_bond].GetOrientation()[i];
  //}
  // sites_[i_bond+1].SubForce(force_);
  // sites_[i_bond].AddForce(force_);
  //}
}

void BeadSpring::UpdateSitePositions(bool midstep) {
  double delta = (midstep ? 0.5 * delta_ : delta_);
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    double const *const f_site = sites_[i_site].GetForce();
    double const *const r_prev = sites_[i_site].GetPrevPosition();
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] = r_prev[i] + f_site[i] * delta / diameter_;
    }
    sites_[i_site].SetPosition(position_);
  }
  // Next, update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site = 0; i_site < n_sites_ - 1; ++i_site) {
    double const *const r_site1 = sites_[i_site].GetPosition();
    double const *const r_site2 = sites_[i_site + 1].GetPosition();
    u_mag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i = 0; i < n_dim_; ++i) r_diff[i] /= u_mag;
    sites_[i_site].SetOrientation(r_diff);
    sites_[i_site].UpdatePeriodic();
  }
  sites_[n_sites_ - 1].SetOrientation(sites_[n_sites_ - 2].GetOrientation());
}

void BeadSpring::UpdateAvgPosition() {
  std::fill(position_, position_ + 3, 0.0);
  std::fill(orientation_, orientation_ + 3, 0.0);
  for (auto site_it = sites_.begin(); site_it != sites_.end(); ++site_it) {
    double const *const site_pos = site_it->GetPosition();
    double const *const site_u = site_it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] += site_pos[i];
      orientation_[i] += site_u[i];
    }
  }
  normalize_vector(orientation_, n_dim_);
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] /= n_sites_;
  }
}

void BeadSpring::ApplyForcesTorques() { ApplyInteractionForces(); }

// FIXME
void BeadSpring::ApplyInteractionForces() {}

// void BeadSpring::Draw(std::vector<graph_struct*> * graph_array) {
// for (auto site=sites_.begin(); site!= sites_.end(); ++site) {
// site->Draw(graph_array);
//}
//}

// Scale bond and site positions from new unit cell
void BeadSpring::ScalePosition() {
  // scale first bond position using new unit cell
  bonds_[0].ScalePosition();
  // then reposition sites based on first bond position
  // handle first site
  double r[3];
  double const *const bond_r = bonds_[0].GetPosition();
  double const *const bond_u = bonds_[0].GetOrientation();
  for (int i = 0; i < n_dim_; ++i)
    r[i] = bond_r[i] - 0.5 * bond_length_ * bond_u[i];
  sites_[0].SetPosition(r);
  // then handle remaining sites
  for (int i_site = 1; i_site < n_sites_; ++i_site) {
    double const *const prev_r = sites_[i_site - 1].GetPosition();
    double const *const prev_u = sites_[i_site - 1].GetOrientation();
    for (int i = 0; i < n_dim_; ++i)
      r[i] = prev_r[i] + bond_length_ * prev_u[i];
    sites_[i_site].SetPosition(r);
  }
  // update remaining bond positions
  UpdateBondPositions();
}

void BeadSpring::GetAvgOrientation(double *au) {
  double avg_u[3] = {0.0, 0.0, 0.0};
  int size = 0;
  for (auto it = sites_.begin(); it != sites_.end(); ++it) {
    double const *const u = it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) avg_u[i] += u[i];
    size++;
  }
  if (size == 0) Logger::Error("Something went wrong in GetAvgOrientation!");
  for (int i = 0; i < n_dim_; ++i) avg_u[i] /= size;
  std::copy(avg_u, avg_u + 3, au);
}

void BeadSpring::GetAvgPosition(double *ap) {
  double avg_p[3] = {0.0, 0.0, 0.0};
  int size = 0;
  for (auto it = sites_.begin(); it != sites_.end(); ++it) {
    double const *const p = it->GetPosition();
    for (int i = 0; i < n_dim_; ++i) avg_p[i] += p[i];
    size++;
  }
  if (size == 0) Logger::Error("Something went wrong in GetAvgPosition!");
  for (int i = 0; i < n_dim_; ++i) avg_p[i] /= size;
  std::copy(avg_p, avg_p + 3, ap);
}

void BeadSpring::ReportAll() {
  printf("cos_thetas:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", cos_thetas_[i]);
  printf("}\n");
}

/* The spec output for one bead_spring is:
    diameter
    length
    persistence_length
    n_sites,
    all site positions
    */

void BeadSpring::WriteSpec(std::fstream &ospec) {
  ospec.write(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  ospec.write(reinterpret_cast<char *>(&length_), sizeof(length_));
  ospec.write(reinterpret_cast<char *>(&persistence_length_),
              sizeof(persistence_length_));
  ospec.write(reinterpret_cast<char *>(&n_sites_), sizeof(n_sites_));
  double temp[3];
  for (int i = 0; i < n_sites_; ++i) {
    double const *const pos = sites_[i].GetPosition();
    std::copy(pos, pos + 3, temp);
    for (auto &p : temp) ospec.write(reinterpret_cast<char *>(&p), sizeof(p));
  }
  return;
}

void BeadSpring::ReadSpec(std::fstream &ispec) {
  if (ispec.eof()) return;
  double r_site[3];
  ispec.read(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  ispec.read(reinterpret_cast<char *>(&length_), sizeof(length_));
  ispec.read(reinterpret_cast<char *>(&persistence_length_),
             sizeof(persistence_length_));
  ispec.read(reinterpret_cast<char *>(&n_sites_), sizeof(n_sites_));
  sites_.resize(n_sites_, sites_[0]);
  // Get initial site position
  // Initialize sites from relative positions
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    for (int i = 0; i < 3; ++i)
      ispec.read(reinterpret_cast<char *>(&r_site[i]), sizeof(double));
    sites_[i_site].SetPosition(r_site);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  UpdatePrevPositions();
  CalculateAngles();
}

void BeadSpring::WritePosit(std::fstream &oposit) {
  double avg_pos[3], avg_u[3];
  GetAvgPosition(avg_pos);
  GetAvgOrientation(avg_u);
  std::copy(avg_pos, avg_pos + 3, position_);
  UpdatePeriodic();
  for (auto &pos : position_)
    oposit.write(reinterpret_cast<char *>(&pos), sizeof(pos));
  for (auto &spos : scaled_position_)
    oposit.write(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : avg_u) oposit.write(reinterpret_cast<char *>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char *>(&length_), sizeof(length_));
}

void BeadSpring::ReadPosit(std::fstream &iposit) {
  if (iposit.eof()) return;
  double avg_pos[3], avg_u[3], s_pos[3];
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&avg_pos[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&s_pos[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&avg_u[i]), sizeof(double));
  iposit.read(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char *>(&length_), sizeof(length_));
  // Initialize first bond position
  for (int i = 0; i < n_dim_; ++i)
    avg_pos[i] = avg_pos[i] - 0.5 * (length_ - bond_length_) * avg_u[i];
  for (int i_bond = 0; i_bond < n_bonds_; ++i_bond) {
    bonds_[i_bond].SetPosition(avg_pos);
    bonds_[i_bond].SetOrientation(avg_u);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].UpdatePeriodic();
    // Set next bond position
    for (int i = 0; i < n_dim_; ++i) avg_pos[i] += bond_length_ * avg_u[i];
  }
}

void BeadSpring::WriteCheckpoint(std::fstream &ocheck) {
  void *rng_state = rng_.GetState();
  size_t rng_size = rng_.GetSize();
  ocheck.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char *>(rng_state), rng_size);
  WriteSpec(ocheck);
}

void BeadSpring::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void *rng_state = rng_.GetState();
  size_t rng_size;
  icheck.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char *>(rng_state), rng_size);
  ReadSpec(icheck);
  n_sites_ = n_bonds_ + 1;
  double r[3];
  for (int i = 0; i < 3; ++i)
    r[i] = bonds_[0].GetPosition()[i] -
           0.5 * bond_length_ * bonds_[0].GetOrientation()[i];
  sites_[0].SetPosition(r);
  sites_[0].SetOrientation(bonds_[0].GetOrientation());
  for (int i_bond = 1; i_bond < n_bonds_; ++i_bond) {
    for (int i = 0; i < 3; ++i)
      r[i] += bond_length_ * sites_[i_bond - 1].GetOrientation()[i];
    sites_[i_bond].SetPosition(r);
    sites_[i_bond].SetOrientation(bonds_[i_bond].GetOrientation());
  }
  for (int i = 0; i < 3; ++i)
    r[i] += bond_length_ * sites_[n_sites_ - 2].GetOrientation()[i];
  sites_[n_sites_ - 1].SetPosition(r);
  sites_[n_sites_ - 1].SetOrientation(bonds_[n_bonds_ - 1].GetOrientation());
  UpdateBondPositions();
  // Reallocate control structures
  cos_thetas_.resize(n_sites_ - 2);  // max_sites-2
}


