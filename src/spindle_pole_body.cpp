#include "spindle_pole_body.h"

// SPB specifics
void SpindlePoleBody::InitConfigurator(const double r,
                                       const double theta,
                                       const double phi,
                                       const double diameter,
                                       const double attach_diameter) {
  if (space_->n_dim != 3) {
    std::cout << "Spindle Pole Bodies require 3 dimensions!\n";
    exit(1);
  }
  diameter_ = diameter;
  attach_diameter_ = attach_diameter;
  std::cout << "Inserting SPB at\n" << std::setprecision(16)
    << "\tr: " << r << ", "
    << "theta: " << theta << ", "
    << "phi: " << phi << "\n";
  conf_rad_ = r;
  double r_anchor[3] = {0.0, 0.0, 0.0};
  r_anchor[0] = r * sin(theta) * cos(phi);
  r_anchor[1] = r * sin(theta) * sin(phi);
  r_anchor[2] = r * cos(theta);
  SetPosition(r_anchor);
  SetPrevPosition(r_anchor);

  UpdateSPBRefVecs();
  SetOrientation(u_anchor_);

  UpdateSPBDragConstants();

  UpdatePeriodic();
}

void SpindlePoleBody::UpdateSPBRefVecs() {
  // Update the u_anchor first (local u pointing to origin)
  for (int i = 0; i < n_dim_; ++i) {
    u_anchor_[i] = -2.0 * position_[i] / conf_rad_;
  }
  double norm_factor = sqrt(1.0/dot_product(3, u_anchor_, u_anchor_));
  for (int i = 0; i < n_dim_; ++i) {
    u_anchor_[i] *= norm_factor;
  }

  double vec0[3] = {1.0, 0.0, 0.0};
  double vec1[3] = {0.0, 1.0, 0.0};
  if (1.0 - ABS(u_anchor_[0]) > 1e-2)
    cross_product(u_anchor_, vec0, v_anchor_, 3);
  else
    cross_product(u_anchor_, vec1, v_anchor_, 3);

  norm_factor = sqrt(1.0/dot_product(3, v_anchor_, v_anchor_));
  for (int i = 0; i < n_dim_; ++i) {
    v_anchor_[i] *= norm_factor;
  }

  cross_product(u_anchor_, v_anchor_, w_anchor_, 3);
}

void SpindlePoleBody::UpdateSPBDragConstants() {
  double spb_diffusion_coefficient = 0.017088; // FIXME

  gamma_tra_ = 1.0 / spb_diffusion_coefficient;
  gamma_rot_ = gamma_tra_*SQR(0.5 * attach_diameter_);
}

void SpindlePoleBody::UpdatePositionMP() {
  // Update the information for the anchor.  There are two torques, the global
  // one maps to translation of the anchor, the local one to the rotation of
  // the anchor (since we're confined to a spherical boundary).  So we have
  // to be careful, as have a translation and a rotation
  //
  // XXX FIXME
  // Most of this is directly taken from bob, so probably can optimize somewhat
  // if needed with correct torques, etc
  double gamma_t = gamma_tra_ * conf_rad_;
  double gamma_r = gamma_rot_;

  double r_anchor_old[3] = {0.0, 0.0, 0.0};
  double u_anchor_old[3] = {0.0, 0.0, 0.0};
  
  double f_dot_u = dot_product(n_dim_, GetForce(), u_anchor_);
  double f_rand[3] = {
    gsl_ran_gaussian_ziggurat(rng_.r, sqrt(2.0 * gamma_t / delta_)),
    gsl_ran_gaussian_ziggurat(rng_.r, sqrt(2.0 * gamma_t / delta_)),
    gsl_ran_gaussian_ziggurat(rng_.r, sqrt(2.0 * gamma_t / delta_))
  };

  // Position
  f_dot_u += dot_product(n_dim_, f_rand, u_anchor_);

  for (int i = 0; i < n_dim_; ++i) {
    u_anchor_old[i] = u_anchor_[i];
    u_anchor_[i] -= (GetForce()[i]+f_rand[i] -
                     f_dot_u * u_anchor_[i]) * delta_ / gamma_t;
  }

  double norm_factor = sqrt(1.0/dot_product(n_dim_, u_anchor_, u_anchor_));
  double r_anchor[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < n_dim_; ++i) {
    u_anchor_[i] *= norm_factor;
    r_anchor_old[i] = position_[i];
    r_anchor[i] = -u_anchor_[i] * conf_rad_;
  }

  SetPrevPosition(r_anchor_old);
  SetPosition(r_anchor);
  AddDr();
  SetOrientation(u_anchor_);

  // Orientation(s)
  double tau_random = gsl_ran_gaussian(rng_.r, sqrt(2.0 * gamma_r / delta_));
  std::copy(torque_, torque_+3, tau_local_);
  double tau_body_full[3] = {
    tau_local_[0] + u_anchor_[0] * tau_random,
    tau_local_[1] + u_anchor_[1] * tau_random,
    tau_local_[2] + u_anchor_[2] * tau_random
  };
  double dv[3] = {0.0, 0.0, 0.0};
  cross_product(tau_body_full, v_anchor_, dv, 3);
  for (int i = 0; i < 3; ++i) dv[i] *= delta_ / gamma_r;
  for (int i = 0; i < 3; ++i) v_anchor_[i] += dv[i];
  double orthocorrection = -dot_product(n_dim_, v_anchor_, u_anchor_);
  for (int i = 0; i < 3; ++i) v_anchor_[i] += orthocorrection * u_anchor_old[i];
  norm_factor = sqrt(1.0/dot_product(n_dim_, v_anchor_, v_anchor_));
  for (int i = 0; i < 3; ++i) v_anchor_[i] *= norm_factor;

  double n[3] = {0.0, 0.0, 0.0};
  cross_product(u_anchor_old, u_anchor_, n, n_dim_);
  double sin_omega = dot_product(n_dim_, n, n);
  norm_factor = 1.0/sqrt(sin_omega);
  if (norm_factor > 1e-8) {
    for (int i = 0; i < 3; ++i) n[i] *= norm_factor;
    double cos_omega = dot_product(n_dim_, u_anchor_, u_anchor_old);
    double vec1[3] = {0.0, 0.0, 0.0};
    cross_product(n, v_anchor_, vec1, n_dim_);
    double factor2 = dot_product(n_dim_, n, v_anchor_) * (1 - cos_omega);
    for (int i = 0; i < 3; ++i) {
      v_anchor_[i] = v_anchor_[i] * cos_omega +
        vec1[i] * sin_omega + n[i] * factor2;
    }
  }

  cross_product(u_anchor_, v_anchor_, w_anchor_, 3);

  // Update the Anchor list!
  UpdateAnchors();
}

void SpindlePoleBody::UpdateAnchors() {
  // find the anchor list
  if (anchors_ == nullptr) return;

  // Get our anchor
  auto malist = (*anchors_)[GetOID()];
  for (auto ait = malist.begin(); ait != malist.end(); ++ait) {
    // Update just the first part, since we are the anchor
    double factor_u = dot_product(n_dim_, ait->pos_rel0_, u_anchor_);
    double factor_v = dot_product(n_dim_, ait->pos_rel0_, v_anchor_);
    double factor_w = dot_product(n_dim_, ait->pos_rel0_, w_anchor_);

    std::cout << "Old Anchor pos0: (" << std::setprecision(16)
      << ait->pos0_[0] << ", " << ait->pos0_[1] << ", " << ait->pos0_[2] << ")\n";
    std::cout << "Old Anchor pos_rel0: (" << std::setprecision(16)
      << ait->pos_rel0_[0] << ", " << ait->pos_rel0_[1] << ", " << ait->pos_rel0_[2] << ")\n";

    // Calculate the new lab frame coordinates
    for (int i = 0; i < 3; ++i) {
      ait->pos0_[i] = position_[i] + factor_u * u_anchor_[i] +
        factor_v * v_anchor_[i] + factor_w * w_anchor_[i];
    }

    double norm_factor = conf_rad_ / sqrt(dot_product(n_dim_, ait->pos0_, ait->pos0_));
    for (int i = 0; i < 3; ++i) ait->pos0_[i] *= norm_factor;

    std::cout << "New Anchor pos0: (" << std::setprecision(16)
      << ait->pos0_[0] << ", " << ait->pos0_[1] << ", " << ait->pos0_[2] << ")\n";
    std::cout << "New Anchor pos_rel0: (" << std::setprecision(16)
      << ait->pos_rel0_[0] << ", " << ait->pos_rel0_[1] << ", " << ait->pos_rel0_[2] << ")\n";
  }
}



// Species specifics
void SpindlePoleBodySpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "SpindlePoleBody species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";
  std::string insertion_type;
  insertion_type = node["spb"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["spb"]["properties"]["overlap"].as<bool>();
  std::cout << "   can overlap:    " << (can_overlap ? "true" : "false") << std::endl;

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 0; // default to orientation
  if (node["spb"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["spb"]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }
  if (node["spb"]["properties"]["draw_type"]) {
    std::string draw_type_s = node["spb"]["properties"]["draw_type"].as<std::string>();
    std::cout << "   draw_type: " << draw_type_s << std::endl;
    if (draw_type_s.compare("flat") == 0) {
      draw_type = 0;
    }
  }

  if (insertion_type.compare("rtp") == 0) {
    if (!can_overlap) {
      std::cout << "Warning, location insertion overrides overlap\n";
      can_overlap = true;
    }
    int nspbs = (int)node["spb"]["spb"].size();
    std::cout << "   nspbs: " << nspbs << std::endl;
    for (int ispb = 0; ispb < nspbs; ++ispb) {
      double diameter         = node["spb"]["spb"][ispb]["diameter"].as<double>();
      double attach_diameter  = node["spb"]["spb"][ispb]["attach_diameter"].as<double>();
      double spb_theta        = node["spb"]["spb"][ispb]["theta"].as<double>();
      double spb_phi          = node["spb"]["spb"][ispb]["phi"].as<double>();
      SpindlePoleBody *member = new SpindlePoleBody(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(0.5 * space_->unit_cell[0][0], spb_theta, spb_phi, diameter, attach_diameter);
      member->SetColor(color, draw_type);
      member->Dump();
      members_.push_back(member);
    }
  } else {
    std::cout << "Insertion type " << insertion_type << " not implemented yet\n";
    exit(1);
  }
}

void SpindlePoleBodySpecies::ConfiguratorSpindle(int ispb, al_set* anchors) {
  char *filename = params_->config_file;
  std::cout << "SpindlePoleBody species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::string insertion_type;
  insertion_type = node["spb"][ispb]["properties"]["insertion_type"].as<std::string>();

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 0; // default to orientation
  if (node["spb"][ispb]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["spb"][ispb]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }
  if (node["spb"][ispb]["properties"]["draw_type"]) {
    std::string draw_type_s = node["spb"][ispb]["properties"]["draw_type"].as<std::string>();
    std::cout << "   draw_type: " << draw_type_s << std::endl;
    if (draw_type_s.compare("flat") == 0) {
      draw_type = 0;
    }
  }

  double diameter         = node["spb"][ispb]["properties"]["diameter"].as<double>();
  double attach_diameter  = node["spb"][ispb]["properties"]["attach_diameter"].as<double>();
  double spb_theta        = node["spb"][ispb]["properties"]["theta"].as<double>();
  double spb_phi          = node["spb"][ispb]["properties"]["phi"].as<double>();

  SpindlePoleBody *member = new SpindlePoleBody(params_, space_, gsl_rng_get(rng_.r), GetSID());
  member->InitConfigurator(0.5 * space_->unit_cell[0][0], spb_theta, spb_phi, diameter, attach_diameter);
  member->SetColor(color, draw_type);
  member->SetAnchors(anchors);
  members_.push_back(member);

  // Create a anchor list entry for this spb
  anchors->insert(std::make_pair(member->GetOID(), std::vector<anchor_t>()));

  member->PrintSPBProperties(ispb);
}
