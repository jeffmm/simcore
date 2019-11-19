#include "rigid_filament.hpp"

RigidFilament::RigidFilament(unsigned long seed) : Mesh(seed) {
  std::fill(body_frame_, body_frame_ + 6, 0.0);
  SetSID(species_id::rigid_filament);
}

void RigidFilament::SetParameters() {
  /* Read parameters from filament parameters */
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  diameter_ = sparams_->diameter;
  max_length_ = sparams_->max_length;
  min_length_ = sparams_->min_length;
  // dynamic_instability_flag_ = sparams_->dynamic_instability_flag;
  // polydispersity_flag_ = sparams_->polydispersity_flag;
  // force_induced_catastrophe_flag_ = sparams_->force_induced_catastrophe_flag;
  // p_g2s_ = sparams_->f_grow_to_shrink * delta_;
  // p_g2p_ = sparams_->f_grow_to_pause * delta_;
  // p_s2p_ = sparams_->f_shrink_to_pause * delta_;
  // p_s2g_ = sparams_->f_shrink_to_grow * delta_;
  // p_p2s_ = sparams_->f_pause_to_shrink * delta_;
  // p_p2g_ = sparams_->f_pause_to_grow * delta_;
  // v_depoly_ = sparams_->v_depoly;
  // v_poly_ = sparams_->v_poly;
  // driving_factor_ = sparams_->driving_factor;
  // stoch_flag_ = params_->stoch_flag;  // include thermal forces
  eq_steps_ = sparams_->n_equil;
  eq_steps_count_ = 0;
  // optical_trap_spring_ = sparams_->optical_trap_spring;
  // optical_trap_flag_ = sparams_->optical_trap_flag;
  // optical_trap_fixed_ = sparams_->optical_trap_fixed;
  // fic_factor_ = sparams_->fic_factor;
  // tip_force_ = 0.0;
  /* Default site in optical trap is the tail */
  // trapped_site_ = 0;

  /* Refine parameters */
}

void RigidFilament::Init(filament_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  InitRigidFilamentLength();
  Reserve();
  InsertRigidFilament(sparams_->insertion_type, -1);
  SetDiffusion();

  // if (optical_trap_flag_) {
  //  double const *const r0 = sites_[0].GetPosition();
  //  std::copy(r0, r0 + 3, optical_trap_pos_);
  //  double const *const r1 = sites_[1].GetPosition();
  //  std::copy(r1, r1 + 3, optical_trap_pos2_);
  //}
  // poly_ = poly_state::grow;
}

/* Returns number of bonds to initialize */
void RigidFilament::InitRigidFilamentLength() {
  if (max_length_ < min_length_) {
    Logger::Warning(
        "Minimum filament length larger than max length -- setting "
        "max_length_ = min_length_");
    max_length_ = min_length_;
  }

  // if (polydispersity_flag_) {
  //  ExponentialDist expon;
  //  expon.Init(length_, min_length_, max_length_);
  //  double roll = rng_.RandomUniform();
  //  length_ = expon.Rand(roll);
  //  if (length_ > max_length_ + 1e-6 || length_ < min_length_ - 1e-6) {
  //    Logger::Error(
  //        "RigidFilament polydispersity distribution returned a value out of "
  //        "range:\nAttempted length: %2.2f, Min length: %2.2f, Max length: "
  //        "%2.2f",
  //        length_, min_length_, max_length_);
  //  }
  //}
  if (length_ > max_length_) {
    Logger::Warning(
        "RigidFilament length larger than max length -- setting length = "
        "max_length");
    length_ = max_length_;
  } else if (length_ < min_length_) {
    Logger::Warning(
        "RigidFilament length less than min length -- setting length = "
        "min_length");
    length_ = min_length_;
  }
  /* Calculate number of bonds in filament at initialization */
  // if (dynamic_instability_flag_) {
  //  n_bonds_ = 2;
  //  bond_length_ = length_ / n_bonds_;
  //  while (dynamic_instability_flag_ && bond_length_ > max_bond_length_) {
  //    n_bonds_ *= 2;
  //    bond_length_ = length_ / n_bonds_;
  //  }
  //} else {
  //  n_bonds_ = (int)floor(length_ / max_bond_length_);
  //  [> Minimum number of bonds must be 2 <]
  //  n_bonds_ = (n_bonds_ == 1 ? 2 : n_bonds_);
  //  bond_length_ = length_ / n_bonds_;
  //}
  //[> Determine maximum number of bonds we may have <]
  // if (dynamic_instability_flag_) {
  //  int max_bonds = (int)ceil(max_length_ / min_bond_length_);
  //  n_bonds_max_ = 2;
  //  while (n_bonds_max_ < max_bonds) {
  //    n_bonds_max_ *= 2;
  //  }
  //} else {
  //  [> Static bond number <]
  //  n_bonds_max_ = n_bonds_;
  //}

  Logger::Trace(
      "RigidFilament initialized with length %2.2f with %d bonds, mesh_id:"
      " %d",
      length_, n_bonds_, GetMeshID());
  true_length_ = length_;
}

void RigidFilament::InsertRigidFilament(std::string insertion_type,
                                        double buffer) {
  if (buffer < 0) {
    buffer = length_;
  }
  if (insertion_type.compare("random") == 0) {
    InsertRandom();
  } else if (insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (insertion_type.compare("centered_random") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    rng_.RandomUnitVector(n_dim_, orientation_);
  } else if (insertion_type.compare("centered_oriented") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else {
    Logger::Error("BrRod insertion type not recognized!");
  }
}

/* Integration scheme taken from Yu-Guo Tao,
   J Chem Phys 122 244903 (2005)
   Explicit calculation of the friction tensor acting on force vector,

   r(t+dt) = r(t) + (Xi^-1 . F_s(t)) * dt + dr(t),

   where friction tensor Xi = gamma_par * |u><u| + gamma_perp * (I - |u><u|),
   u is the orientation, F_s is the force from interactions, and dr(t) is the
   random displacement due to random force,

   dr(t) = Xi^-1 . F_r * dt,

   which is treated separately as a random displacement with std dev
   sqrt(2*kT*dt/gamma_(par/perp)) along par/perp unit vectors
   relative to rod. */
void RigidFilament::Integrate() {
  // Explicit calculation of Xi.F_s
  double dr[3] = {0, 0, 0};
  for (int i = 0; i < n_dim_; ++i) {
    for (int j = 0; j < n_dim_; ++j) {
      dr[i] +=
          gamma_par_ * orientation_[i] * orientation_[j] * force_[j] * delta_;
    }
    dr[i] += force_[i] * gamma_perp_ * delta_;
    position_[i] += dr[i];
  }
  // Reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, 3);  // ndim=3 since torques
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] += du[i] * delta_ / gamma_rot_;
  }
  // Add the random displacement dr(t)
  AddRandomDisplacement();
  // Update the orientation due to torques and random rotation
  AddRandomReorientation();
  UpdatePeriodic();
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void RigidFilament::AddRandomDisplacement() {
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  double mag = rng_.RandomNormal(diffusion_par_);
  for (int i = 0; i < n_dim_; ++i) position_[i] += mag * orientation_[i];
  // Then the perpendicular component(s)
  for (int j = 0; j < n_dim_ - 1; ++j) {
    mag = rng_.RandomNormal(diffusion_perp_);
    for (int i = 0; i < n_dim_; ++i)
      position_[i] += mag * body_frame_[n_dim_ * j + i];
  }
  // Handle the random orientation update after updating orientation from
  // interaction torques
}

/* The orientation update is also from Yu-Guo Tao,

   u(t+dt) = u(t) + gamma_rot^-1 * T_s(t) x u(t) * dt + du(t)

   where similar to above, du(t) is the reorientation due to
   random forces, and is treated as random displacement vector(s)
   orthogonal to u(t) with std dev sqrt(2*kT*dt/gamma_rot) */
void RigidFilament::AddRandomReorientation() {
  // Now handle the random orientation update
  for (int j = 0; j < n_dim_ - 1; ++j) {
    double mag = rng_.RandomNormal(diffusion_rot_);
    for (int i = 0; i < n_dim_; ++i) {
      orientation_[i] += mag * body_frame_[n_dim_ * j + i];
    }
  }
  normalize_vector(orientation_, n_dim_);
}

void RigidFilament::GetBodyFrame() {
  if (n_dim_ == 2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  } else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]), n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

void RigidFilament::InsertAt(const double *const new_pos,
                             const double *const u) {
  Logger::Trace(
      "Inserting filament at [%2.1f, %2.1f, %2.1f] with orientation"
      "[%2.1f, %2.1f, %2.1f]",
      new_pos[0], new_pos[1], new_pos[2], u[0], u[1], u[2]);
  RelocateMesh(new_pos, u);
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
  // if (optical_trap_flag_) {
  //  trapped_site_ = 0;
  /* For cilia flag, if filament is oriented along a negative
   * dimension, assume the filament needs to be fixed at the plus
   * end */
  //  if (cilia_trap_flag_) {
  //    for (int i = 0; i < n_dim_; ++i) {
  //      if (u[i] < 0) {
  //        trapped_site_ = n_sites_ - 1;
  //      }
  //    }
  //  }
  //  int trapped_2 = (trapped_site_ == 0 ? 1 : n_sites_ - 2);
  //  const double *const r0 = sites_[trapped_site_].GetPosition();
  //  const double *const r1 = sites_[trapped_2].GetPosition();
  //  std::copy(r0, r0 + 3, optical_trap_pos_);
  //  std::copy(r1, r1 + 3, optical_trap_pos2_);
  //}
  // poly_ = poly_state::grow;
}

void RigidFilament::SetDiffusion() {
  // Sets the diffusion in accordance to Lowen, Phys. Rev. E, 1994.
  double L = length_ + diameter_;
  double p = L / diameter_;
  double log_p = log(p);
  gamma_par_ = 2.0 * L / 3.0 / (log_p - 0.207 + 0.980 / p - 0.133 / SQR(p));
  gamma_perp_ = 4.0 * L / 3.0 / (log_p + 0.839 + 0.185 / p + 0.233 / SQR(p));
  gamma_rot_ = CUBE(L) / 9.0 / (log_p - 0.662 + 0.917 / p - 0.050 / SQR(p));
  diffusion_par_ = sqrt(2 * delta_ / gamma_par_);
  diffusion_perp_ = sqrt(2 * delta_ / gamma_perp_);
  diffusion_rot_ = sqrt(2 * delta_ / gamma_rot_);
}

double const RigidFilament::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_ * length_ + 0.25 * M_PI * diameter_ * diameter_;
  } else {
    return 0.25 * M_PI * diameter_ * diameter_ * length_ +
           1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_;
  }
}

void RigidFilament::UpdatePosition() {
  ApplyForcesTorques();
  Integrate();
  UpdateAvgPosition();
  // DynamicInstability();
  eq_steps_count_++;
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
// void RigidFilament::Integrate() {
//  CalculateAngles();
//  CalculateTangents();
//  if (midstep_) {
//    ConstructUnprojectedRandomForces();
//    GeometricallyProjectRandomForces();
//    UpdatePrevPositions();
//  }
//  AddRandomForces();
//  UpdateSitePositions();
//  UpdateBondPositions();
//}

void RigidFilament::CalculateTangents() {
  for (auto it = sites_.begin(); it != sites_.end(); ++it) {
    it->CalcTangent();
  }
}

void RigidFilament::GetNematicOrder(double *nematic_order_tensor) {
  for (auto it = bonds_.begin(); it != bonds_.end(); ++it) {
    double const *const u_bond = it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        if (i == j) {
          nematic_order_tensor[3 * i + j] +=
              (2 * u_bond[i] * u_bond[j] - 1) / n_bonds_;
        } else {
          nematic_order_tensor[3 * i + j] +=
              2 * u_bond[i] * u_bond[j] / n_bonds_;
        }
      }
    }
  }
}

void RigidFilament::GetPolarOrder(double *polar_order_vector) {
  double p_i[3] = {0, 0, 0};
  for (auto it = bonds_.begin(); it != bonds_.end(); ++it) {
    double const *const u_bond = it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      p_i[i] += u_bond[i];
    }
  }
  for (int i = 0; i < n_dim_; ++i) {
    polar_order_vector[i] += p_i[i] / n_bonds_;
  }
}

void RigidFilament::AddRandomForces() {
  if (!stoch_flag_) return;
  for (auto site = sites_.begin(); site != sites_.end(); ++site)
    site->AddRandomForce();
}

void RigidFilament::UpdateSitePositions() {
  double f_site[3];
  // First get total forces
  // Handle end sites first
  // Now update positions
  double f_term[3], r_new[3];
  int site_index = 0;
  int next_site = n_dim_ * n_dim_;
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
  }
  sites_[n_sites_ - 1].SetOrientation(sites_[n_sites_ - 2].GetOrientation());
  // Finally, normalize site positions, making sure the sites are still
  // rod-length apart
  // if (CheckBondLengths()) {
  //  for (int i_site = 1; i_site < n_sites_; ++i_site) {
  //    double const *const r_site1 = sites_[i_site - 1].GetPosition();
  //    double const *const u_site1 = sites_[i_site - 1].GetOrientation();
  //    for (int i = 0; i < n_dim_; ++i)
  //      r_diff[i] = r_site1[i] + bond_length_ * u_site1[i];
  //    sites_[i_site].SetPosition(r_diff);
  //  }
  //}
}

bool RigidFilament::CheckBondLengths() {
  bool renormalize = false;
  for (int i_site = 1; i_site < n_sites_; ++i_site) {
    double const *const r_site1 = sites_[i_site - 1].GetPosition();
    double const *const r_site2 = sites_[i_site].GetPosition();
    double a = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      double temp = r_site2[i] - r_site1[i];
      a += temp * temp;
    }
    a = sqrt(a);
    double err = ABS(bond_length_ - a) / bond_length_;
    if (err > 1e-3) {
      renormalize = true;
    }
  }
  return renormalize;
}

void RigidFilament::UpdateAvgPosition() {
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

void RigidFilament::ApplyForcesTorques() {
  ApplyInteractionForces();
  // if (optical_trap_flag_) {
  //  double f_trap1[3] = {0};
  //  double f_trap2[3] = {0};
  //  double const *const r0 = sites_[trapped_site_].GetPosition();
  //  int trap2 = (trapped_site_ == 0 ? 1 : n_sites_ - 2);
  //  double const *const r1 = sites_[trap2].GetPosition();
  //  for (int i = 0; i < n_dim_; ++i) {
  //    f_trap1[i] = optical_trap_spring_ * (optical_trap_pos_[i] - r0[i]);
  //    if (optical_trap_fixed_) {
  //      f_trap2[i] = optical_trap_spring_ * (optical_trap_pos2_[i] - r1[i]);
  //    }
  //  }
  //  sites_[trapped_site_].AddForce(f_trap1);
  //  if (optical_trap_fixed_) {
  //    sites_[trap2].AddForce(f_trap2);
  //  }
  //}
  // if (anchored_) ApplyAnchorForces();
}

void RigidFilament::ApplyInteractionForces() {
  double pure_torque[3] = {0, 0, 0};
  double site_force[3] = {0, 0, 0};
  double linv = 1.0 / bond_length_;
  if (!sparams_->drive_from_bond_center) {
    // Driving originating from the site tangents
    CalculateTangents();
  }
  for (int i = 0; i < n_bonds_; ++i) {
    double const *const f = bonds_[i].GetForce();
    double const *const t = bonds_[i].GetTorque();
    double const *const u = sites_[i].GetOrientation();
    // if (i == n_bonds_ - 1) {
    //  tip_force_ = -dot_product(n_dim_, u, f);
    //}
    AddPotential(bonds_[i].GetPotentialEnergy());
    // Convert torques into forces at bond ends
    // u x t / bond_length = pure torque force at tail of bond
    cross_product(u, t, pure_torque, 3);
    for (int i = 0; i < n_dim_; ++i) {
      pure_torque[i] *= linv;
      site_force[i] = 0.5 * f[i];
    }
    // Add translational forces and pure torque forces at bond ends
    sites_[i].AddForce(site_force);
    sites_[i].AddForce(pure_torque);
    for (int j = 0; j < n_dim_; ++j) pure_torque[j] *= -1;
    sites_[i + 1].AddForce(site_force);
    sites_[i + 1].AddForce(pure_torque);
    // The driving factor is a force per unit length,
    // so need to multiply by bond length to get f_dr on bond
    // if (eq_steps_count_ > eq_steps_) {
    //  double f_dr[3] = {};
    //  if (sparams_->drive_from_bond_center) {
    //    // Add driving (originating from the com of the bond)
    //    double mag = 0.5 * driving_factor_ * bond_length_;
    //    for (int j = 0; j < n_dim_; ++j) f_dr[j] = mag * u[j];
    //    sites_[i].AddForce(f_dr);
    //    sites_[i + 1].AddForce(f_dr);
    //  } else {
    //    // Driving from sites
    //    double mag = length_ * driving_factor_ / n_sites_;
    //    double const *const u_tan = sites_[i].GetTangent();
    //    for (int j = 0; j < n_dim_; ++j) {
    //      f_dr[j] = mag * u_tan[j];
    //    }
    //    sites_[i].AddForce(f_dr);
    //  }
    //}
  }
}

// void RigidFilament::DynamicInstability() {
//  if (midstep_ || !dynamic_instability_flag_) return;
//  UpdatePolyState();
//  GrowRigidFilament();
//  SetDiffusion();
//}

// void RigidFilament::GrowRigidFilament() {
//  // If the filament is paused, do nothing
//  if (poly_ == +poly_state::pause) return;
//  // Otherwise, adjust filament length due to polymerization
//  double delta_length = 0;
//  if (poly_ == +poly_state::grow) {
//    delta_length = v_poly_ * delta_;
//  } else if (poly_ == +poly_state::shrink) {
//    if (n_bonds_ == 2 && bond_length_ <= min_bond_length_) {
//      return;
//    }
//    delta_length = -v_depoly_ * delta_;
//  }
//  length_ += delta_length;
//  RescaleBonds();
//  if (bond_length_ > max_bond_length_) {
//    DoubleGranularityLinear();
//    // RebindMotors();
//  } else if (bond_length_ < min_bond_length_ && n_bonds_ > 2) {
//    HalfGranularityLinear();
//    // RebindMotors();
//  }
//}

// void RigidFilament::UpdatePolyState() {
//  double p_g2s = p_g2s_;
//  double p_p2s = p_p2s_;
//  double roll = rng_.RandomUniform();
//  double p_norm;
//  // Modify catastrophe probabilities if the end of the filament is under a
//  // load
//  if (force_induced_catastrophe_flag_ && tip_force_ > 0) {
//    double p_factor = exp(fic_factor_ * tip_force_);
//    p_g2s = (p_g2s + p_g2p_) * p_factor;
//    p_p2s = p_p2s * p_factor;
//  }
//  // RigidFilament shrinking
//  if (poly_ == +poly_state::shrink) {
//    p_norm = p_s2g_ + p_s2p_;
//    if (p_norm > 1.0)
//      poly_ = (roll < p_s2g_ / p_norm ? poly_state::grow : poly_state::pause);
//    else {
//      if (roll < p_s2g_)
//        poly_ = poly_state::grow;
//      else if (roll < (p_s2g_ + p_s2p_))
//        poly_ = poly_state::pause;
//    }
//  }
//  // RigidFilament growing
//  else if (poly_ == +poly_state::grow) {
//    p_norm = p_g2s + p_g2p_;
//    if (p_norm > 1.0)
//      poly_ = (roll < p_g2s / p_norm ? poly_state::shrink :
//      poly_state::pause);
//    else {
//      if (roll < p_g2s)
//        poly_ = poly_state::shrink;
//      else if (roll < (p_g2s + p_g2p_))
//        poly_ = poly_state::pause;
//    }
//  }
//  // RigidFilament paused
//  else if (poly_ == +poly_state::pause) {
//    p_norm = p_p2g_ + p_p2s;
//    if (p_norm > 1)
//      poly_ = (roll < p_p2g_ / p_norm ? poly_state::grow :
//      poly_state::shrink);
//    else {
//      if (roll < p_p2g_)
//        poly_ = poly_state::grow;
//      else if (roll < (p_p2g_ + p_p2s))
//        poly_ = poly_state::shrink;
//    }
//  }
//  // Check to make sure the filament lengths stay in the correct ranges
//  if (length_ < min_length_)
//    poly_ = poly_state::grow;
//  else if (length_ > max_length_)
//    poly_ = poly_state::shrink;
//}

// void RigidFilament::CheckFlocking() {
//  double avg_polar_order = 0;
//  double avg_contact_number = 0;
//  for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
//    [> Check if filament satisfies the flock condition by checking whether
//       the average bond polar order is above cutoff. Also see whether the
//       average contact number exceeds the cutoff for interior/exterior
//       filament. */
//    avg_polar_order += bond->GetPolarOrder();
//    avg_contact_number += bond->GetContactNumber();
//  }
//  avg_polar_order /= n_bonds_;
//  avg_contact_number /= n_bonds_;

//  int in_flock_prev = in_flock_;
//  in_flock_ = 0;
//  flock_change_state_ = 0;
//  if (avg_polar_order >= params_->flock_polar_min) {
//    // RigidFilament is in a flock
//    if (in_flock_prev == 0) {
//      // RigidFilament joined flock this timestep
//      flock_change_state_ = 1;
//    }
//    if (avg_contact_number >= params_->flock_contact_min) {
//      // RigidFilament is in flock interior
//      in_flock_ = 1;
//    } else {
//      // RigidFilament is in flock exterior
//      in_flock_ = 2;
//    }
//  } else if (in_flock_prev > 0) {
//    // RigidFilament left flock this timestep
//    flock_change_state_ = 2;
//  }
//}

void RigidFilament::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    bond->SetFlockType(in_flock_);
    bond->Draw(graph_array);
  }
}

// Scale bond and site positions from new unit cell
void RigidFilament::ScalePosition() {
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

void RigidFilament::ReportAll() {
  // printf("tensions:\n  {");
  // for (int i = 0; i < n_sites_ - 1; ++i) printf(" %5.5f ", tensions_[i]);
  // printf("}\n");
  // printf("cos_thetas:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", cos_thetas_[i]);
  // printf("}\n");
  // printf("g_mat_lower:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", g_mat_lower_[i]);
  // printf("}\n");
  // printf("g_mat_upper:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", g_mat_upper_[i]);
  // printf("}\n");
  // printf("g_mat_diag:\n  {");
  // for (int i = 0; i < n_sites_ - 1; ++i) printf(" %5.5f ", g_mat_diag_[i]);
  // printf("}\n");
  // printf("det_t_mat:\n  {");
  // for (int i = 0; i < n_sites_ + 1; ++i) printf(" %5.5f ", det_t_mat_[i]);
  // printf("}\n");
  // printf("det_b_mat:\n  {");
  // for (int i = 0; i < n_sites_ + 1; ++i) printf(" %5.5f ", det_b_mat_[i]);
  // printf("}\n");
  // printf("h_mat_diag:\n  {");
  // for (int i = 0; i < n_sites_ - 1; ++i) printf(" %5.5f ", h_mat_diag_[i]);
  // printf("}\n");
  // printf("h_mat_upper:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", h_mat_upper_[i]);
  // printf("}\n");
  // printf("h_mat_lower:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", h_mat_lower_[i]);
  // printf("}\n");
  // printf("k_eff:\n  {");
  // for (int i = 0; i < n_sites_ - 2; ++i) printf(" %5.5f ", k_eff_[i]);
  // printf("}\n\n\n");
}

/* The spec output for one filament is:
    diameter
    length
    persistence_length (added 1/17/2017)
    friction_par (added 1/17/2017)
    friction_perp (added 1/17/2017)
    bond_length
    n_bonds,
    position of first site
    position of last site
    all bond orientations
    */

void RigidFilament::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing rigid filament specs, object id: %d", GetOID());
  Mesh::WriteSpec(ospec);
  // ospec.write(reinterpret_cast<char *>(&persistence_length_),
  // sizeof(double)); ospec.write(reinterpret_cast<char *>(&poly_),
  // sizeof(unsigned char));
}

/* double diameter
   double length
   double persistence length
   double friction_par
   double friction_perp
   double bond_length
   double n_bonds
   double[3] pos_tail
   double[3] pos_head
   double[n_bonds*3] bond_orientations
*/
void RigidFilament::ReadSpec(std::fstream &ispec) {
  if (ispec.eof()) return;
  Mesh::ReadSpec(ispec);
  // ispec.read(reinterpret_cast<char *>(&persistence_length_), sizeof(double));
  // ispec.read(reinterpret_cast<char *>(&poly_), sizeof(unsigned char));
  // CalculateAngles();
}

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void RigidFilament::WritePosit(std::fstream &oposit) {
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

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void RigidFilament::ReadPosit(std::fstream &iposit) {
  if (iposit.eof()) return;
  posits_only_ = true;
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
  std::copy(avg_pos, avg_pos + 3, position_);
  for (int i = 0; i < n_dim_; ++i)
    avg_pos[i] = avg_pos[i] - 0.5 * length_ * avg_u[i];
  for (int i_bond = 0; i_bond < n_bonds_; ++i_bond) {
    sites_[i_bond].SetPosition(avg_pos);
    sites_[i_bond].SetOrientation(avg_u);
    for (int i = 0; i < n_dim_; ++i) {
      avg_pos[i] += 0.5 * bond_length_ * avg_u[i];
    }
    bonds_[i_bond].SetPosition(avg_pos);
    bonds_[i_bond].SetOrientation(avg_u);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].UpdatePeriodic();
    // Set next bond position
    for (int i = 0; i < n_dim_; ++i) {
      avg_pos[i] += 0.5 * bond_length_ * avg_u[i];
    }
  }
  sites_[n_bonds_].SetPosition(avg_pos);
  sites_[n_bonds_].SetOrientation(avg_u);
  SetOrientation(avg_u);
  UpdatePeriodic();
  CalculateAngles();
}

void RigidFilament::WriteCheckpoint(std::fstream &ocheck) {
  Mesh::WriteCheckpoint(ocheck);
}

void RigidFilament::ReadCheckpoint(std::fstream &icheck) {
  Mesh::ReadCheckpoint(icheck);
}
