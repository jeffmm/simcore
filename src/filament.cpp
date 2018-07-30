#include "filament.h"

Filament::Filament() : Mesh() {
  SetParameters();
} 

void Filament::SetParameters() {
  color_ = params_->filament.color;
  draw_ = draw_type::_from_string(params_->filament.draw_type.c_str());
  length_ = params_->filament.length;
  persistence_length_ = params_->filament.persistence_length;
  if (params_->filament.perlen_ratio > 0) {
    persistence_length_ = length_ * params_->filament.perlen_ratio;
  }
  diameter_ = params_->filament.diameter;
  // TODO JMM: add subdivisions of bonds for interactions, 
  //           should depend on cell length
  max_length_ = params_->filament.max_length;
  min_length_ = params_->filament.min_length;
  max_bond_length_ = params_->filament.max_bond_length;
  min_bond_length_ = params_->filament.min_bond_length;
  if (length_ > 0 && params_->filament.n_bonds > 0) {
    min_bond_length_ = length_ / params_->filament.n_bonds + 1e-6;
    max_bond_length_ = length_ / params_->filament.n_bonds - 1e-6;
  }
  dynamic_instability_flag_ = params_->filament.dynamic_instability_flag;
  spiral_flag_ = params_->filament.spiral_flag;
  force_induced_catastrophe_flag_ = params_->filament.force_induced_catastrophe_flag;
  p_g2s_ = params_->filament.f_grow_to_shrink*delta_;
  p_g2p_ = params_->filament.f_grow_to_pause*delta_;
  p_s2p_ = params_->filament.f_shrink_to_pause*delta_;
  p_s2g_ = params_->filament.f_shrink_to_grow*delta_;
  p_p2s_ = params_->filament.f_pause_to_shrink*delta_;
  p_p2g_ = params_->filament.f_pause_to_grow*delta_;
  v_depoly_ = params_->filament.v_depoly;
  v_poly_ = params_->filament.v_poly;
  driving_factor_ = params_->filament.driving_factor;
  friction_ratio_ = params_->filament.friction_ratio;
  metric_forces_ = params_->filament.metric_forces;
  stoch_flag_ = params_->stoch_flag; // determines whether we are using stochastic forces
  eq_steps_ = params_->n_steps_equil;
  eq_steps_count_ = 0;
  fic_factor_ = params_->filament.fic_factor;
  tip_force_ = 0.0;
  k_off_ = params_->motor.k_off;
  k_on_ = params_->motor.k_on;
  motor_concentration_ = params_->motor.concentration;
  n_motors_bound_ = 0;
  shuffle_flag_ = params_->filament.shuffle;
  shuffle_factor_ = params_->filament.shuffle_factor;
  shuffle_frequency_ = params_->filament.shuffle_frequency;
  if (params_->filament.spiral_number_fail_condition <= 0) {
    params_->filament.spiral_number_fail_condition = 1e-6;
  }
}

void Filament::SetAnchor(Anchor * a) {
  anchor_ = a;
  anchored_ = true;
}

bool Filament::CheckBounds(double buffer) {
  for (auto site=sites_.begin();site!=sites_.end(); ++site) {
    if (site->CheckBounds()) return true;
  }
  return false;
}

void Filament::Init(bool force_overlap) {
  InitElements();
  InsertFilament(force_overlap);
  if (spiral_flag_) {
    InitSpiral2D();
    return;
  }
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
  poly_ = poly_state::grow;
}

void Filament::InitElements() {
  int n_bonds;
  if (max_length_ < min_length_) {
    warning("Minimum filament length larger than max length -- setting min_length = max_length");
    min_length_ = max_length_;
  }
  if (length_ > max_length_) {
    warning("Filament length larger than max length -- setting length = max_length");
    length_ = max_length_;
  }
  else if (length_ > 0 && length_ < min_length_) {
    warning("Filament length less than min length -- setting length = min_length");
    length_ = min_length_;
  }
  int max_bonds = (int) ceil(max_length_/min_bond_length_);
  n_bonds_max_ = 2;
  while (n_bonds_max_ < max_bonds) {
    n_bonds_max_*=2;
  }

  //if (max_length_/n_bonds_max_ < min_bond_length_) {
    //error_exit("min_length_ of flexible filament bonds too large for filament length.");
  //}
  // Initialize mesh
  Reserve(n_bonds_max_);
  //Allocate control structures
  int n_sites_max = n_bonds_max_+1;
  tensions_.resize(n_sites_max-1); //max_sites -1
  g_mat_lower_.resize(n_sites_max-2); //max_sites-2
  g_mat_upper_.resize(n_sites_max-2); //max_sites-2
  g_mat_diag_.resize(n_sites_max-1); //max_sites-1
  det_t_mat_.resize(n_sites_max+1); //max_sites+1
  det_b_mat_.resize(n_sites_max+1); //max_sites+1
  g_mat_inverse_.resize(n_sites_max-2); //max_sites-2
  k_eff_.resize(n_sites_max-2); //max_sites-2
  h_mat_diag_.resize(n_sites_max-1); //max_sites-1
  h_mat_upper_.resize(n_sites_max-2); //max_sites-2
  h_mat_lower_.resize(n_sites_max-2); //max_sites-2
  gamma_inverse_.resize(n_sites_max*n_dim_*n_dim_); //max_sites*ndim*ndim
  cos_thetas_.resize(n_sites_max-2); //max_sites-2
}

bool Filament::InsertFirstBond() {
  if (anchored_) {
    InitSiteAt(anchor_->position_,diameter_);
    std::copy(anchor_->orientation_,anchor_->orientation_+3,orientation_);
    GenerateProbableOrientation();
    AddBondToTip(orientation_,bond_length_);
  }
  else if (params_->filament.insertion_type.compare("random") == 0) {
    InitRandomSite(diameter_);
    AddRandomBondToTip(bond_length_);
  }
  else if (params_->filament.insertion_type.compare("random_nematic") == 0) {
    InitRandomSite(diameter_);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = (gsl_rng_uniform_pos(rng_.r) > 0.5 ? 1.0 : -1.0);
    AddBondToTip(orientation_, bond_length_);
  }
  else if (params_->filament.insertion_type.compare("random_polar") == 0) {
    InitRandomSite(diameter_);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
    AddBondToTip(orientation_, bond_length_);
  }
  else if (params_->filament.insertion_type.compare("centered_oriented") == 0) {
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
    for (int i=0;i<n_dim_;++i) {
      position_[i] = - 0.5*length_*orientation_[i];
    }
    InitSiteAt(position_,diameter_);
    AddBondToTip(orientation_, bond_length_);
  }
  else if (params_->filament.insertion_type.compare("centered_random") == 0) {
    generate_random_unit_vector(n_dim_,orientation_,rng_.r);
    for (int i=0;i<n_dim_;++i) {
      position_[i] = - 0.5*length_*orientation_[i];
    }
    InitSiteAt(position_,diameter_);
    AddBondToTip(orientation_, bond_length_);
  }
  else {
    // Assume custom arrangement for now
    InitSiteAt(position_,diameter_);
    AddBondToTip(orientation_, bond_length_);
    return false;
  }
  return true;
}

void Filament::InsertFilament(bool force_overlap) {
  //if (params_->filament.insertion_type.compare("custom") == 0) {
    //return;
  //}
  bool out_of_bounds = true;
  int n_bonds;
  bool polydisperse = (length_<=0);
  do {
    out_of_bounds = false;
    // Polydisperse filaments
    if (polydisperse) {
      length_ = min_length_ + (max_length_-min_length_)*gsl_rng_uniform_pos(rng_.r);
    }
    //do {
    n_bonds = floor(length_/max_bond_length_);
    if (dynamic_instability_flag_) {
      n_bonds = 2;
    }
    bond_length_ = length_/n_bonds;
    if (dynamic_instability_flag_) {
      while (bond_length_ > max_bond_length_) {
        n_bonds *= 2;
        bond_length_ = length_/n_bonds;
      }
    }
    //if (n_bonds == 2 && bond_length_ <= diameter_) {
      //error_exit("bond_length <= diameter despite minimum number of bonds.\nTry reducing filament diameter or increasing filament length.");
    //}
      //if (bond_length_ <= diameter_) {
        //max_bond_length_ += 0.1*max_bond_length_;
        //warning("bond_length <= diameter, increasing max_bond_length to %2.2f",max_bond_length_);
      //}
    //} while (bond_length_ <= diameter_);
    Clear();
    force_overlap = !InsertFirstBond();
    if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) continue;
    SetOrientation(bonds_[n_bonds_-1].GetOrientation());
    bool probable_orientation = (params_->filament.insertion_type.compare("simple_crystal") != 0 && params_->filament.insertion_type.compare("random_polar") != 0 && params_->filament.insertion_type.compare("random_nematic") != 0);
    for (int i=0;i<n_bonds-1;++i) {
      if (probable_orientation) {
        GenerateProbableOrientation();
      }
      AddBondToTip(orientation_, bond_length_);
      if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) break;
    }
  } while (out_of_bounds);
  for (bond_iterator bond=bonds_.begin(); bond!=bonds_.end(); ++bond) {
    bond->SetColor(color_,draw_);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  shuffle_factor_ = 2*(gsl_rng_uniform_pos(rng_.r)-0.5)*shuffle_factor_;
}

void Filament::InsertAt(double *pos, double *u) {
  std::copy(pos,pos+3,position_);
  std::copy(u,u+3,orientation_);
  normalize_vector(orientation_,n_dim_);
  for (int i=0;i<n_dim_;++i) {
    position_[i] -= 0.5*length_*orientation_[i];
  }

  bool force_overlap = true;
  bool out_of_bounds = true;
  int n_bonds;
  bool polydisperse = (length_<=0);
  do {
    out_of_bounds = false;
    // Polydisperse filaments
    if (polydisperse) {
      length_ = min_length_ + (max_length_-min_length_)*gsl_rng_uniform_pos(rng_.r);
    }
    //do {
    n_bonds = floor(length_/max_bond_length_);
    if (dynamic_instability_flag_) {
      n_bonds = 2;
    }
    bond_length_ = length_/n_bonds;
    if (dynamic_instability_flag_) {
      while (bond_length_ > max_bond_length_) {
        n_bonds *= 2;
        bond_length_ = length_/n_bonds;
      }
    }
    Clear();
    InitSiteAt(position_,diameter_);
    AddBondToTip(orientation_, bond_length_);
    //InsertFirstBond();
    if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) continue;
    SetOrientation(bonds_[n_bonds_-1].GetOrientation());
    for (int i=0;i<n_bonds-1;++i) {
      AddBondToTip(orientation_, bond_length_);
      if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) break;
    }
  } while (out_of_bounds);
  for (bond_iterator bond=bonds_.begin(); bond!=bonds_.end(); ++bond) {
    bond->SetColor(color_,draw_);
  }
  shuffle_factor_ = 2*(gsl_rng_uniform_pos(rng_.r)-0.5)*shuffle_factor_;
  UpdateBondPositions();
  UpdateSiteOrientations();
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
  poly_ = poly_state::grow;
}

// Place a spool centered at the origin
void Filament::InitSpiral2D() {
  if (n_dim_ > 2)
    error_exit("3D Spirals not coded yet.");
  double prev_pos[3] = {0, 0, 0};
  std::fill(position_, position_+3, 0);
  sites_[n_sites_-1].SetPosition(prev_pos);
  for (auto site=sites_.begin(); site!= sites_.end(); ++site) {
    site->SetDiameter(diameter_);
    site->SetLength(bond_length_);
  }
  double step = diameter_ / M_PI;
  double theta = bond_length_ / step;
  for (int i=2; i<n_sites_+1; ++i) {
    double move = step*theta;
    double angle = theta + 2*M_PI;
    position_[0] = move * cos(angle);
    position_[1] = move * sin(angle);
    theta += bond_length_ / move;
    // Set current site position and orientation
    sites_[n_sites_-i].SetPosition(position_);
    for (int j=0; j<3; ++j) {
      orientation_[j] = (prev_pos[j] - position_[j])/bond_length_;
      prev_pos[j] = position_[j];
    }
    sites_[n_sites_-i].SetOrientation(orientation_);
  }
  // Set last site orientation
  sites_[n_sites_-1].SetOrientation(sites_[n_sites_-2].GetOrientation());
  // Init usual suspects
  UpdatePrevPositions();
  CalculateAngles();
  UpdateBondPositions();
  UpdateSiteOrientations();
  SetDiffusion();
  poly_ = poly_state::grow;
}

void Filament::SetDiffusion() {
  double logLD = log(length_/diameter_);
  //double gamma_0 = 4.0/3.0*eps*((1+0.64*eps)/(1-1.15*eps) + 1.659 * SQR(eps));
  //friction_perp_ = bond_length_ * gamma_0;
  friction_perp_ = 4.0*length_/(3.0*n_sites_*logLD);
  friction_par_ = friction_perp_ / friction_ratio_;
  rand_sigma_perp_ = sqrt(24.0*friction_perp_/delta_);
  rand_sigma_par_ = sqrt(24.0*friction_par_/delta_);
}

void Filament::GenerateProbableOrientation() {
  /* This updates the current orientation with a generated probable orientation where
  we generate random theta pulled from probability distribution P(th) = exp(k cos(th))
  where k is the persistence_length of the filament
  If k is too large, there is enormous imprecision in this calculation since sinh(k)
  is very large so to fix this I introduce an approximate distribution that is valid 
  for large k */
  double theta;
  if (persistence_length_ == 0) {
    theta = gsl_rng_uniform_pos(rng_.r) * M_PI;
  }
  else if (persistence_length_ < 100) {
    theta = gsl_rng_uniform_pos(rng_.r) * M_PI;
    theta = acos( log( exp(-persistence_length_/bond_length_) + 
          2.0*gsl_rng_uniform_pos(rng_.r)*sinh(persistence_length_/bond_length_) ) 
        / (persistence_length_/bond_length_) );
  }
  else {
    theta = acos( (log( 2.0*gsl_rng_uniform_pos(rng_.r)) - 
          log(2.0) + persistence_length_/bond_length_)
        /(persistence_length_/bond_length_) );
  }
  double new_orientation[3] = {0, 0, 0};
  if (n_dim_==2) {
    theta = (gsl_rng_uniform_int(rng_.r,2)==0 ? -1 : 1) * theta;
    new_orientation[0] = cos(theta);
    new_orientation[1] = sin(theta);
  }
  else {
    double phi = gsl_rng_uniform_pos(rng_.r) * 2.0 * M_PI;
    new_orientation[0] = sin(theta)*cos(phi);
    new_orientation[1] = sin(theta)*sin(phi);
    new_orientation[2] = cos(theta);
  }
  rotate_vector_relative(n_dim_, new_orientation, orientation_);
  std::copy(new_orientation,new_orientation+3,orientation_);
}

double const Filament::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_*length_ + 0.25*M_PI*diameter_*diameter_;
  }
  if (n_dim_ == 3) {
    return 0.25*M_PI*diameter_*diameter_*length_ + 1.0/6.0 * M_PI * diameter_*diameter_*diameter_;
  }
}

void Filament::UpdatePosition(bool midstep) {
  midstep_ = midstep;
  if (!midstep_) {
    CalculateBinding();
    for (auto it=motors_.begin(); it!= motors_.end(); ++it) {
      it->UpdatePosition();
    }
  }
  ApplyForcesTorques();
  Integrate();
  // FIXME temporary
  UpdateAvgPosition();
  DynamicInstability();
  eq_steps_count_++;
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
void Filament::Integrate() {
  CalculateAngles();
  CalculateSpiralNumber();
  CalculateTangents();
  if (midstep_) {
    ConstructUnprojectedRandomForces();
    GeometricallyProjectRandomForces();
    UpdatePrevPositions();
  }
  AddRandomForces();
  CalculateBendingForces();
  CalculateTensions();
  UpdateSitePositions();
  UpdateBondPositions();
  UpdateSiteOrientations();
}

//  FIXME move this out of filament!
void Filament::CalculateBinding() {
  if (++n_step_ % 100 != 0) return;
  double vol = GetVolume();
  // Check motor unbinding
  if (gsl_rng_uniform_pos(rng_.r) <= k_off_*n_motors_bound_*100*delta_) {
    UnbindMotor();
  }
  // Check motor binding
  if (gsl_rng_uniform_pos(rng_.r) <= motor_concentration_*vol*k_on_*100*delta_) {
    BindMotor();
  }
}

void Filament::BindMotor() {
  Motor m;
  motors_.push_back(m);
  motors_.back().Init();
  int i_bond = gsl_rng_uniform_int(rng_.r,n_bonds_);
  motors_.back().AttachBondRandom(&bonds_[i_bond],i_bond*bond_length_);
  n_motors_bound_++;
  motor_concentration_ -= 1.0/space_->volume;
}

void Filament::UnbindMotor() {
  int i_motor = gsl_rng_uniform_int(rng_.r,n_motors_bound_);
  motors_.erase(motors_.begin()+i_motor);
  n_motors_bound_--;
  motor_concentration_ += 1.0/space_->volume;
}

void Filament::UpdateSiteOrientations() {
  for (int i=0; i<n_sites_-1; ++i) {
    sites_[i].SetOrientation(bonds_[i].GetOrientation());
  }
  sites_[n_sites_-1].SetOrientation(bonds_[n_bonds_-1].GetOrientation());
}

void Filament::CalculateAngles(bool rescale) {
  double cos_angle, angle_sum = 0;
  bool sharp_angle = false;
  for (int i_site=0; i_site<n_sites_-2; ++i_site) {
    double const * const u1 = sites_[i_site].GetOrientation();
    double const * const u2 = sites_[i_site+1].GetOrientation();
    cos_angle = dot_product(n_dim_, u1, u2);
    cos_thetas_[i_site] = cos_angle;
    if (spiral_flag_) {
      angle_sum += acos(cos_angle);
    }
  }
  if (rescale && sharp_angle && midstep_ && length_/(n_bonds_+1) > min_bond_length_) {
    DoubleGranularityLinear();
    RebindMotors();
    UpdateSiteOrientations();
    CalculateAngles(false);
  }
}

void Filament::CalculateTangents() {
  for (auto it = sites_.begin(); it!=sites_.end(); ++it) {
    it->CalcTangent();
  }
}

void Filament::CalculateSpiralNumber() {
  // This will calculate the angle 
  double u_bond[3] = {0,0,0};
  double const * const u_head = bonds_[n_bonds_-1].GetOrientation();
  double const * const p_head = sites_[n_sites_-1].GetPosition();
  //double const * const p_head_real = sites_[n_sites_-1].GetPosition();
   //From COM
  //double p_head[3] = {0,0,0};
  //GetAvgPosition(p_head);
  //for (int i=0;i<n_dim_;++i) {
    //u_bond[i] = p_head_real[i] - p_head[i];
  //}
  //normalize_vector(u_bond,n_dim_);
  std::copy(u_head,u_head+3,u_bond);
  for (int i=0;i<3;++i) {
    u_bond[i] = -u_bond[i];
  }
  spiral_number_ = 0;
  for (int i=n_sites_-3; i>=0; --i) {
  //for (int i=n_sites_-2; i>=0; --i) {
    double u[3] = {0,0,0};
    double const * const p = sites_[i].GetPosition();
    for (int i=0;i<n_dim_;++i) {
      u[i] = p[i] - p_head[i];
    }
    normalize_vector(u,n_dim_);
    double dp = dot_product(n_dim_,u_bond,u);
    /* Note that floating point errors occasionally lead to 
       NANs in the acos function, these fix that */
    if (dp > 1) dp = 1;
    else if (dp < -1) dp = -1;
    double angle = acos(dp);
    double temp[3];
    cross_product(u_bond,u,temp,3);
    int sign = SIGNOF(temp[2]);
    for (int i=0;i<n_dim_; ++i) {
      u_bond[i] = u[i];
    }
    spiral_number_ += sign*angle;
  }
  // Failed spiral (waiting until the filament straightens)
  if (spiral_flag_ && ABS(GetSpiralNumber()) < params_->filament.spiral_number_fail_condition) {
    early_exit = true;
  }
}

double Filament::GetSpiralNumber() {
  return 0.5*spiral_number_/M_PI;
  //double temp = 0.5*spiral_number_/M_PI;
  //if (temp >= 0.75) return temp+0.25;
  //else return temp/0.75;
  //return temp;
}

void Filament::GetNematicOrder(double * nematic_order_tensor) {
  for (auto it=bonds_.begin(); it!=bonds_.end(); ++it) {
    double const * const u_bond = it->GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      for (int j=0; j<n_dim_; ++j) {
        if (i==j) {
          nematic_order_tensor[3*i+j] += (2*u_bond[i]*u_bond[j] - 1)/n_bonds_;
        }
        else {
          nematic_order_tensor[3*i+j] += 2*u_bond[i]*u_bond[j]/n_bonds_;
        }
      }
    }
  }
}

void Filament::GetPolarOrder(double * polar_order_vector) {
  double p_i[3] = {0,0,0};
  for (auto it=bonds_.begin(); it!=bonds_.end(); ++it) {
    double const * const u_bond = it->GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      p_i[i] += u_bond[i];
    }
  }
  for (int i=0; i<n_dim_; ++i) {
    polar_order_vector[i] += p_i[i]/n_bonds_;
  }
}

void Filament::ConstructUnprojectedRandomForces() {
  // Create unprojected forces, see J. Chem. Phys. 122, 084903 (2005), eqn. 40.
  // xi is the random force vector with elements that are uncorrelated and randomly
  // distributed uniformly between -0.5 and 0.5, xi_term is the outer product of the
  // tangent vector u_tan_i u_tan_i acting on the vector xi
  if (!stoch_flag_) return;
  double xi[3], xi_term[3], f_rand[3];
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    double const * const utan = sites_[i_site].GetTangent();
    for (int i=0; i<n_dim_; ++i) 
      xi[i] = gsl_rng_uniform_pos(rng_.r) - 0.5;
    if (n_dim_ == 2) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0];
    }
    else if (n_dim_ == 3) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1]
                  + utan[0] * utan[2] * xi[2];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0]
                  + utan[1] * utan[2] * xi[2];
      xi_term[2] = SQR(utan[2]) * xi[2] + utan[0] * utan[2] * xi[0]
                  + utan[1] * utan[2] * xi[1];
    }
    for (int i=0; i<n_dim_; ++i) {
      f_rand[i] = rand_sigma_perp_ * xi[i] + 
          (rand_sigma_par_ - rand_sigma_perp_) * xi_term[i];
    }
    sites_[i_site].SetRandomForce(f_rand);
  }
}

void Filament::GeometricallyProjectRandomForces() {
  if (!stoch_flag_) return;
  double f_rand_temp[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    // Use the tensions vector to calculate the hard components of the random forces
    // These are not the same as the tensions, they will be calculated later
    double const * const f_rand1 = sites_[i_site].GetRandomForce();
    double const * const f_rand2 = sites_[i_site+1].GetRandomForce();
    double const * const u_site = sites_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      f_rand_temp[i] = f_rand2[i] - f_rand1[i];
    tensions_[i_site] = dot_product(n_dim_, f_rand_temp, u_site);
    // Then get the G arrays (for inertialess case where m=1, see
    // ref. 15 of above paper)
    g_mat_diag_[i_site] = 2;
    if (i_site > 0) {
      g_mat_upper_[i_site-1] = - cos_thetas_[i_site-1];
      g_mat_lower_[i_site-1] = - cos_thetas_[i_site-1];
    }
  }
  // Now solve using tridiagonal solver
  tridiagonal_solver(&g_mat_lower_, &g_mat_diag_, &g_mat_upper_, &tensions_, n_sites_-1);
  // Update to the projected brownian forces
  // First the end sites:
  double f_proj[3];
  for (int i=0; i<n_dim_; ++i) {
    f_proj[i] = sites_[0].GetRandomForce()[i]
      + tensions_[0] * sites_[0].GetOrientation()[i];
  }
  sites_[0].SetRandomForce(f_proj);
  for (int i=0; i<n_dim_; ++i) {
    f_proj[i] = sites_[n_sites_-1].GetRandomForce()[i] 
    - tensions_[n_sites_-2] * sites_[n_sites_-2].GetOrientation()[i];
  }
  sites_[n_sites_-1].SetRandomForce(f_proj);
  // Then the rest
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    double const * const u1 = sites_[i_site-1].GetOrientation();
    double const * const u2 = sites_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      f_proj[i] = sites_[i_site].GetRandomForce()[i] + tensions_[i_site]*u2[i] - tensions_[i_site-1] * u1[i];
    }
    sites_[i_site].SetRandomForce(f_proj);
  }
}

void Filament::AddRandomForces() {
  if (!stoch_flag_) return;
  for (auto site=sites_.begin(); site!=sites_.end(); ++site)
    site->AddRandomForce();
}

void Filament::CalculateBendingForces() {
  if (metric_forces_) {
    det_t_mat_[0] = 1;
    det_t_mat_[1] = 2;
    det_b_mat_[n_sites_] = 1;
    det_b_mat_[n_sites_-1] = 2;
    for (int i=2; i<n_sites_; ++i) {
      det_t_mat_[i] = 2 * det_t_mat_[i-1] - SQR(- cos_thetas_[i-2]) * det_t_mat_[i-2];
      det_b_mat_[n_sites_-i] = 2 * det_b_mat_[n_sites_-i+1] - SQR(- cos_thetas_[n_sites_-i-1]) * det_b_mat_[n_sites_-i+2];
    }
    double det_g = det_t_mat_[n_sites_-1];
    for(int i=0; i<n_sites_-2; ++i) {
      g_mat_inverse_[i] = cos_thetas_[i] * det_t_mat_[i] * det_b_mat_[i+3] / det_g;
    }
  }
  else {
    for(int i=0; i<n_sites_-2; ++i) {
      g_mat_inverse_[i] = 0;
    }
  }                                        
  // Now calculate the effective rigidities 
  for (int i=0; i<n_sites_-2; ++i) {
    k_eff_[i] = (persistence_length_ + bond_length_ * g_mat_inverse_[i])/SQR(bond_length_);
  }
  // Using these, we can calculate the forces on each of the sites
  // These calculations were done by hand and are not particularly readable,
  // but are more efficient than doing it explicitly in the code for readability
  // If this ever needs fixed, you need to either check the indices very carefully
  // or redo the calculation by hand! 
  // See Pasquali and Morse, J. Chem. Phys. Vol 116, No 5 (2002)
  double f_site[3] = {0, 0, 0};
  if (n_dim_ == 2) {
    for (int k_site=0; k_site<n_sites_; ++k_site) {
      std::fill(f_site,f_site+3,0.0);
      if (k_site>1) {
        double const * const u1 = sites_[k_site-2].GetOrientation();
        double const * const u2 = sites_[k_site-1].GetOrientation();
        f_site[0] += k_eff_[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] );
        f_site[1] += k_eff_[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        double const * const u1 = sites_[k_site-1].GetOrientation();
        double const * const u2 = sites_[k_site].GetOrientation();
        f_site[0] += k_eff_[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1]
                        -((1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1]) );
        f_site[1] += k_eff_[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0]
                        -((1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0]) );
      }
      if (k_site<n_sites_-2) {
        double const * const u1 = sites_[k_site].GetOrientation();
        double const * const u2 = sites_[k_site+1].GetOrientation();
        f_site[0] -= k_eff_[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] );
        f_site[1] -= k_eff_[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0] );
      }
      sites_[k_site].AddForce(f_site);
    }
  }
  else if (n_dim_ == 3) {
    for (int k_site=0; k_site<n_sites_; ++k_site) {
      std::fill(f_site,f_site+3,0.0);
      if (k_site>1) {
        double const * const u1 = sites_[k_site-2].GetOrientation();
        double const * const u2 = sites_[k_site-1].GetOrientation();
        f_site[0] += k_eff_[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] );
        f_site[1] += k_eff_[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] );
        f_site[2] += k_eff_[k_site-2] * ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[2]*u2[1]*u1[1] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        double const * const u1 = sites_[k_site-1].GetOrientation();
        double const * const u2 = sites_[k_site].GetOrientation();
        f_site[0] += k_eff_[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] 
                        - ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] ) );
        f_site[1] += k_eff_[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2]
                        - ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] ) );
        f_site[2] += k_eff_[k_site-1] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1]
                        - ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[1]*u2[2]*u1[1] ) );
      }
      if(k_site<n_sites_-2) {
        double const * const u1 = sites_[k_site].GetOrientation();
        double const * const u2 = sites_[k_site+1].GetOrientation();
        f_site[0] -= k_eff_[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] );
        f_site[1] -= k_eff_[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2] );
        f_site[2] -= k_eff_[k_site] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1] );
      }
      sites_[k_site].AddForce(f_site);
    }
  }
}

void Filament::CalculateTensions() {
  // Calculate friction_inverse matrix
  int site_index = 0;
  int next_site = n_dim_*n_dim_;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    int gamma_index = 0;
    double const * const utan = sites_[i_site].GetTangent();
    for (int i=0; i<n_dim_; ++i) {
      for (int j=0; j<n_dim_; ++j) {
        gamma_inverse_[site_index+gamma_index] = 1.0/friction_par_ * (utan[i]*utan[j])
                                     + 1.0/friction_perp_ * ((i==j ? 1 : 0) - utan[i]*utan[j]);
        gamma_index++;
      }
    }
    site_index += next_site;
  }
  // Populate the H matrix and Q vector using tensions (p_vec) array
  double temp_a, temp_b;
  double f_diff[3];
  double utan1_dot_u2, utan2_dot_u2;
  site_index = 0;
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    // f_diff is the term in par_entheses in equation 29 of J. Chem. Phys. 122, 084903 (2005)
    double const * const f1 = sites_[i_site].GetForce();
    double const * const f2 = sites_[i_site+1].GetForce();
    double const * const u2 = sites_[i_site].GetOrientation();
    double const * const utan1 = sites_[i_site].GetTangent();
    double const * const utan2 = sites_[i_site+1].GetTangent();
    for (int i=0; i<n_dim_; ++i) {
      temp_a = gamma_inverse_[site_index+n_dim_*i] * f1[0] + gamma_inverse_[site_index + n_dim_*i+1] * f1[1];
      if (n_dim_ == 3) temp_a += gamma_inverse_[site_index+n_dim_*i+2] * f1[2];
      temp_b = gamma_inverse_[site_index+next_site+n_dim_*i] * f2[0] + gamma_inverse_[site_index+next_site+n_dim_*i+1] * f2[1];
      if (n_dim_ == 3) temp_b += gamma_inverse_[site_index+next_site+n_dim_*i+2] * f2[2];
      f_diff[i] = temp_b - temp_a;
    }
    tensions_[i_site] = dot_product(n_dim_, u2, f_diff);
    utan1_dot_u2 = dot_product(n_dim_, utan1, u2);
    utan2_dot_u2 = dot_product(n_dim_, utan2, u2);
    h_mat_diag_[i_site] = 2.0/friction_perp_ + (1.0/friction_par_ - 1.0/friction_perp_) *
      (SQR(utan1_dot_u2) + SQR(utan2_dot_u2));
    if (i_site>0) {
      double const * const u1 = sites_[i_site-1].GetOrientation();
      h_mat_upper_[i_site-1] = -1.0/friction_perp_ * dot_product(n_dim_, u2, u1)
        - (1.0/friction_par_ - 1.0/friction_perp_) * (dot_product(n_dim_, utan1, u1) *
            dot_product(n_dim_, utan1, u2));
      h_mat_lower_[i_site-1] = h_mat_upper_[i_site-1];
    }
    site_index += next_site;
  }
  tridiagonal_solver(&h_mat_lower_, &h_mat_diag_, &h_mat_upper_, &tensions_, n_sites_-1);
}

void Filament::UpdateSitePositions() {
  double delta = (midstep_ ? 0.5*delta_ : delta_);
  double f_site[3];
  // First get total forces
  // Handle end sites first
  for (int i=0; i<n_dim_; ++i)
    f_site[i] = tensions_[0] * sites_[0].GetOrientation()[i];
  sites_[0].AddForce(f_site);
  for (int i=0; i<n_dim_; ++i)
    f_site[i] = -tensions_[n_sites_-2] * sites_[n_sites_-2].GetOrientation()[i];
  sites_[n_sites_-1].AddForce(f_site);
  // and then the rest
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    double const * const u_site1 = sites_[i_site-1].GetOrientation();
    double const * const u_site2 = sites_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      f_site[i] = tensions_[i_site] * u_site2[i] - tensions_[i_site-1] * u_site1[i];
    }
    sites_[i_site].AddForce(f_site);
  }
  // Now update positions
  double f_term[3], r_new[3];
  int site_index = 0;
  int next_site = n_dim_*n_dim_;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    double const * const f_site1 = sites_[i_site].GetForce();
    double const * const r_site1 = sites_[i_site].GetPosition();
    double const * const r_prev = sites_[i_site].GetPrevPosition();
    for (int i=0; i<n_dim_; ++i) {
      f_term[i] = gamma_inverse_[site_index+n_dim_*i] * f_site1[0] + gamma_inverse_[site_index+n_dim_*i+1] * f_site1[1];
      if (n_dim_ == 3)
        f_term[i] += gamma_inverse_[site_index+n_dim_*i+2] * f_site1[2];
      r_new[i] = r_prev[i] + delta * f_term[i];
    }
    sites_[i_site].SetPosition(r_new);
    site_index += next_site;
  }
  // Next, update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    double const * const r_site1 = sites_[i_site].GetPosition();
    double const * const r_site2 = sites_[i_site+1].GetPosition();
    u_mag = 0.0;
    for (int i=0; i<n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i=0; i<n_dim_; ++i)
      r_diff[i]/=u_mag;
    sites_[i_site].SetOrientation(r_diff);
  }
  sites_[n_sites_-1].SetOrientation(sites_[n_sites_-2].GetOrientation());
  // Finally, normalize site positions, making sure the sites are still rod-length apart
  if (CheckBondLengths()) {
    for (int i_site=1; i_site<n_sites_; ++i_site) {
      double const * const r_site1 = sites_[i_site-1].GetPosition();
      double const * const u_site1 = sites_[i_site-1].GetOrientation();
      for (int i=0; i<n_dim_; ++i)
        r_diff[i] = r_site1[i] + bond_length_ * u_site1[i];
      sites_[i_site].SetPosition(r_diff);
    }
  }
}

bool Filament::CheckBondLengths() {
  bool renormalize = false;
  for (int i_site=1;i_site<n_sites_;++i_site) {
    double const * const r_site1 = sites_[i_site-1].GetPosition();
    double const * const r_site2 = sites_[i_site].GetPosition();
    double a = 0.0;
    for (int i=0;i<n_dim_;++i) {
      double temp = r_site2[i] - r_site1[i];
      a += temp*temp;
    }
    a = sqrt(a);
    double err = ABS(bond_length_ - a)/bond_length_;
    if (err > 1e-3) {
      renormalize = true;
    }
  }
  return renormalize;
}

void Filament::UpdateAvgPosition() {
  std::fill(position_, position_+3, 0.0);
  std::fill(orientation_, orientation_+3, 0.0);
  for (auto site_it=sites_.begin(); site_it!=sites_.end(); ++site_it) {
    double const * const site_pos = site_it->GetPosition();
    double const * const site_u = site_it->GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      position_[i] += site_pos[i];
      orientation_[i] += site_u[i];
    }
  }
  normalize_vector(orientation_, n_dim_);
  for (int i=0; i<n_dim_; ++i) {
    position_[i] /= n_sites_;
  }
}

void Filament::ApplyForcesTorques() {
  ApplyInteractionForces();
  if (anchored_) ApplyAnchorForces();
}

void Filament::ApplyInteractionForces() {
  double pure_torque[3] = {0,0,0};
  double site_force[3] = {0,0,0};
  double linv=1.0/bond_length_;
  for (int i=0; i<n_bonds_; ++i) {
    double const * const f = bonds_[i].GetForce();
    double const * const t = bonds_[i].GetTorque();
    double const * const u = sites_[i].GetOrientation();
    if (i == n_bonds_-1) {
      tip_force_ = -dot_product(n_dim_,u,f);
    }
    AddPotential(bonds_[i].GetPotentialEnergy());
    // Convert torques into forces at bond ends
    // u x t / bond_length = pure torque force at tail of bond
    cross_product(u, t, pure_torque, 3);
    for (int i=0; i<n_dim_; ++i) {
      pure_torque[i]*=linv;
      site_force[i] = 0.5*f[i];
    }
    // Add translational forces and pure torque forces at bond ends
    sites_[i].AddForce(site_force);
    sites_[i].AddForce(pure_torque);
    for (int j=0; j<n_dim_; ++j)
      pure_torque[j] *= -1;
    sites_[i+1].AddForce(site_force);
    sites_[i+1].AddForce(pure_torque);
    // The driving factor is a force per unit length,
    // so need to multiply by bond length to get f_dr on bond
    if (eq_steps_count_ > eq_steps_) {
      if (params_->filament.driving_method == 0) {
        // Add driving (originating from the com of the bond)
        double f_dr[3];
        for (int j=0; j<n_dim_; ++j)
          f_dr[j] = 0.5*u[j]*driving_factor_ * bond_length_;
        sites_[i].AddForce(f_dr);
        sites_[i+1].AddForce(f_dr);
      }
    }
    else if (shuffle_flag_) {
      if (gsl_rng_uniform_pos(rng_.r) < shuffle_frequency_) {
        shuffle_factor_ = 2*(gsl_rng_uniform_pos(rng_.r)-0.5)*params_->filament.shuffle_factor;
      }
      double f_dr[3];
      for (int j=0; j<n_dim_; ++j) {
        f_dr[j] = 0.5*u[j]*shuffle_factor_ * bond_length_;
      }
      sites_[i].AddForce(f_dr);
      sites_[i+1].AddForce(f_dr);
    }
  }
  if (params_->filament.driving_method == 1) {
    // Add driving (originating from the site tangents)
    CalculateTangents();
    double mag = length_/n_sites_*driving_factor_;
    for (int i=0; i<n_sites_; ++i) {
      double f_dr[3];
      double const * const u_tan = sites_[i].GetTangent();
      for (int j=0; j<n_dim_; ++j) {
        f_dr[j] = mag*u_tan[j];
      }
      sites_[i].AddForce(f_dr);
    }
  }
}

void Filament::ApplyAnchorForces() {
  double const * const tail_pos = sites_[0].GetPosition();
  double dr[3] = {0,0,0};
  double dr_mag=0.0;
  for (int i=0;i<n_dim_;++i) {
    dr[i] = tail_pos[i] - anchor_->position_[i];
    dr_mag += SQR(dr[i]);
  }
  dr_mag = sqrt(dr_mag);
  if (dr_mag < 1e-12) {
    std::fill(anchor_->force_,anchor_->force_+3,0.0);
  }
  else {
    for (int i=0;i<n_dim_;++i) {
      anchor_->force_[i] = anchor_->k_spring_*(dr_mag - anchor_->spring_length_)*dr[i]/dr_mag;
    }
  }
  sites_[0].SubForce(anchor_->force_);

  if (anchor_->alignment_potential_) {
    double const * const tail_u = bonds_[0].GetOrientation();
    double cos_theta = dot_product(n_dim_,tail_u,anchor_->orientation_);
    double factor = anchor_->k_align_* sqrt(persistence_length_)*0.5*ABS(1-cos_theta);
    double temp[3] = {0,0,0};
    cross_product(anchor_->orientation_, tail_u, temp, 3);
    normalize_vector(temp, 3);
    for (int i=0;i<3;++i) {
      anchor_->torque_[i] = factor * temp[i]; 
    }
    bonds_[0].SubTorque(anchor_->torque_);
  }
}

void Filament::DynamicInstability() {
  if (midstep_ || !dynamic_instability_flag_) return;
  UpdatePolyState();
  GrowFilament();
  SetDiffusion();
}

void Filament::GrowFilament() {
  // If the filament is paused, do nothing
  if (poly_ == +poly_state::pause) return;
  // Otherwise, adjust filament length due to polymerization
  double delta_length = 0;
  if (poly_ == +poly_state::grow) {
    delta_length = v_poly_ * delta_;
  }
  else if (poly_ == +poly_state::shrink) {
    delta_length = -v_depoly_ * delta_;
  }
  length_ += delta_length;
  RescaleBonds();
  if (bond_length_ > max_bond_length_) {
    DoubleGranularityLinear();
    RebindMotors();
    UpdateSiteOrientations();
  }
  else if (bond_length_ < min_bond_length_ && n_bonds_ > 2) {
    HalfGranularityLinear();
    RebindMotors();
    UpdateSiteOrientations();
  }
}

void Filament::RescaleBonds() {
  UpdatePrevPositions();
  double old_bond_length = bond_length_;
  bond_length_ = length_/n_bonds_;
  double k, dl;
  double r2[3] = {0,0,0};
  if (poly_ == +poly_state::shrink) {
    // Old code
    double const * const r0=GetTailPosition();
    double const * const u0=GetTailOrientation();
    for (int i=0; i<n_dim_; ++i) {
      r2[i] = r0[i] + u0[i] * bond_length_;
    }
    sites_[1].SetPosition(r2);
    dl = old_bond_length - bond_length_;
    for (int i_site=2; i_site<n_sites_; ++i_site) {
      double const * const r1 = sites_[i_site-1].GetPrevPosition();
      double const * const u1 = sites_[i_site-1].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-2]) - SQR(dl) + SQR(bond_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site-2]*dl + sqrt(k);
      for (int i=0; i<n_dim_; ++i)
        r2[i] = r1[i] + k * u1[i];
      sites_[i_site].SetPosition(r2);
      dl = old_bond_length - k;
    }
  }
  else if (poly_ == +poly_state::grow) {
    dl = old_bond_length;
    for (int i_site=1; i_site<n_sites_-1; ++i_site) {
      double const * const r_old = sites_[i_site].GetPrevPosition();
      double const * const u = sites_[i_site].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-1]) - SQR(dl) + SQR(bond_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site-1]*dl + sqrt(k);
      for (int i=0; i<n_dim_; ++i) {
        r2[i] = r_old[i] + k * u[i];
      }
      sites_[i_site].SetPosition(r2);
      dl = old_bond_length - k;
    }
    double const * const u = sites_[n_sites_-2].GetOrientation();
    double const * const r_old = sites_[n_sites_-2].GetPosition();
    for (int i=0; i<n_dim_; ++i) {
      r2[i] = r_old[i] + bond_length_ * u[i];
    }
    sites_[n_sites_-1].SetPosition(r2);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  CalculateAngles(false);
}


void Filament::UpdatePolyState() {
  double p_g2s = p_g2s_;
  double p_p2s = p_p2s_;
  double roll = gsl_rng_uniform_pos(rng_.r);
  double p_norm;
  // Modify catastrophe probabilities if the end of the filament is under a load
  if (force_induced_catastrophe_flag_ && tip_force_ > 0) {
    double p_factor = exp(fic_factor_*tip_force_);
    p_g2s = (p_g2s+p_g2p_)*p_factor;
    p_p2s = p_p2s*p_factor;
  }
  // Filament shrinking
  if (poly_ == +poly_state::shrink) {
    p_norm = p_s2g_ + p_s2p_;
    if (p_norm > 1.0) 
      poly_ = (roll < p_s2g_/p_norm ? poly_state::grow : poly_state::pause);
    else {
      if (roll < p_s2g_) 
        poly_ = poly_state::grow;
      else if (roll < (p_s2g_ + p_s2p_)) 
        poly_ = poly_state::pause;
    }
  }
  // Filament growing
  else if (poly_ == +poly_state::grow) {
    p_norm = p_g2s + p_g2p_;
    if (p_norm > 1.0)
      poly_ = (roll < p_g2s/p_norm ? poly_state::shrink : poly_state::pause);
    else {
      if (roll < p_g2s) 
        poly_ = poly_state::shrink;
      else if (roll < (p_g2s + p_g2p_)) 
        poly_ = poly_state::pause;
    }
  }
  // Filament paused
  else if (poly_ == +poly_state::pause) {
    p_norm = p_p2g_ + p_p2s;
    if (p_norm > 1) 
      poly_ = (roll < p_p2g_/p_norm ? poly_state::grow : poly_state::shrink);
    else {
      if (roll < p_p2g_) 
        poly_ = poly_state::grow;
      else if (roll < (p_p2g_ + p_p2s)) 
        poly_ = poly_state::shrink;
    }
  }
  // Check to make sure the filament lengths stay in the correct ranges
  if (length_ < min_length_)
    poly_ = poly_state::grow;
  else if (length_ > max_length_)
    poly_ = poly_state::shrink;
}

void Filament::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto bond=bonds_.begin(); bond!= bonds_.end(); ++bond) {
    bond->Draw(graph_array);
  }
  // FIXME temporary
  for (auto motor=motors_.begin(); motor!= motors_.end(); ++motor) {
    motor->Draw(graph_array);
  }
}

// Scale bond and site positions from new unit cell
void Filament::ScalePosition() {
  // scale first bond position using new unit cell
  bonds_[0].ScalePosition();
  // then reposition sites based on first bond position
  // handle first site
  double r[3];
  double const * const bond_r = bonds_[0].GetPosition();
  double const * const bond_u = bonds_[0].GetOrientation();
  for (int i=0; i<n_dim_; ++i)
    r[i] = bond_r[i] - 0.5*bond_length_*bond_u[i];
  sites_[0].SetPosition(r);
  // then handle remaining sites
  for (int i_site=1; i_site<n_sites_; ++i_site) {
    double const * const prev_r = sites_[i_site-1].GetPosition();
    double const * const prev_u = sites_[i_site-1].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      r[i] = prev_r[i] + bond_length_*prev_u[i];
    sites_[i_site].SetPosition(r);
  }
  // update remaining bond positions
  UpdateBondPositions();
}


void Filament::GetAvgOrientation(double * au) {
  double avg_u[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=sites_.begin(); it!=sites_.end(); ++it) {
    double const * const u = it->GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      avg_u[i] += u[i];
    size++;
  }
  if (size == 0)
    error_exit("Something went wrong in GetAvgOrientation!");
  for (int i=0; i<n_dim_; ++i)
    avg_u[i]/=size;
  std::copy(avg_u, avg_u+3, au);
}

void Filament::RebindMotors() {
  for (motor_iterator it=motors_.begin(); it!=motors_.end(); ++it) {
    double mesh_lambda = it->GetMeshLambda();
    int i_site = (int) floor(mesh_lambda/bond_length_);
    double bond_lambda = mesh_lambda - i_site*bond_length_;
    if (i_site > n_sites_-2) {
      i_site = n_sites_-2;
      mesh_lambda = length_;
      bond_lambda = bond_length_;
    }
    if (i_site < 0) {
      i_site = 0;
      mesh_lambda = 0;
      bond_lambda = 0;
    }
    it->AttachToBond(sites_[i_site].GetOutgoingBond(),bond_lambda,mesh_lambda);
  }
}

void Filament::GetAvgPosition(double * ap) {
  double avg_p[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=sites_.begin(); it!=sites_.end(); ++it) {
    double const * const p = it->GetPosition();
    for (int i=0; i<n_dim_; ++i)
      avg_p[i] += p[i];
    size++;
  }
  if (size == 0)
    error_exit("Something went wrong in GetAvgPosition!");
  for (int i=0; i<n_dim_; ++i)
    avg_p[i]/=size;
  std::copy(avg_p, avg_p+3, ap);
}

void Filament::ReportAll() {
  printf("tensions:\n  {");
  for (int i=0; i<n_sites_-1; ++i)
    printf(" %5.5f ",tensions_[i]);
  printf("}\n");
  printf("cos_thetas:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",cos_thetas_[i]);
  printf("}\n");
  printf("g_mat_lower:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",g_mat_lower_[i]);
  printf("}\n");
  printf("g_mat_upper:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",g_mat_upper_[i]);
  printf("}\n");
  printf("g_mat_diag:\n  {");
  for (int i=0; i<n_sites_-1; ++i)
    printf(" %5.5f ",g_mat_diag_[i]);
  printf("}\n");
  printf("det_t_mat:\n  {");
  for (int i=0; i<n_sites_+1; ++i)
    printf(" %5.5f ",det_t_mat_[i]);
  printf("}\n");
  printf("det_b_mat:\n  {");
  for (int i=0; i<n_sites_+1; ++i)
    printf(" %5.5f ",det_b_mat_[i]);
  printf("}\n");
  printf("h_mat_diag:\n  {");
  for (int i=0; i<n_sites_-1; ++i)
    printf(" %5.5f ",h_mat_diag_[i]);
  printf("}\n");
  printf("h_mat_upper:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",h_mat_upper_[i]);
  printf("}\n");
  printf("h_mat_lower:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",h_mat_lower_[i]);
  printf("}\n");
  printf("k_eff:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",k_eff_[i]);
  printf("}\n\n\n");
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

void Filament::WriteSpec(std::fstream &ospec){
  ospec.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ospec.write(reinterpret_cast<char*>(&length_), sizeof(length_));
  ospec.write(reinterpret_cast<char*>(&persistence_length_), sizeof(persistence_length_));
  ospec.write(reinterpret_cast<char*>(&friction_par_), sizeof(friction_par_));
  ospec.write(reinterpret_cast<char*>(&friction_perp_), sizeof(friction_perp_));
  ospec.write(reinterpret_cast<char*>(&bond_length_), sizeof(bond_length_));
  ospec.write(reinterpret_cast<char*>(&n_bonds_), sizeof(n_bonds_));
  double temp[3];
  double const * const r0 = sites_[0].GetPosition();
  std::copy(r0, r0+3, temp);
  for (auto& pos : temp)
    ospec.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  double const * const rf = sites_[n_bonds_].GetPosition();
  std::copy(rf, rf+3, temp);
  for (auto& pos : temp)
    ospec.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for (int i=0; i<n_bonds_; ++i) {
    double const * const orientation = bonds_[i].GetOrientation();
    std::copy(orientation, orientation+3, temp);
    for (auto& u : temp) 
      ospec.write(reinterpret_cast<char*>(&u), sizeof(u));
  }
  return;
}

void Filament::ReadSpec(std::fstream &ispec) {
  if (ispec.eof()) return;
  double r0[3], rf[3], u_bond[3];
  ispec.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ispec.read(reinterpret_cast<char*>(&length_), sizeof(length_));
  ispec.read(reinterpret_cast<char*>(&persistence_length_), sizeof(persistence_length_));
  ispec.read(reinterpret_cast<char*>(&friction_par_), sizeof(friction_par_));
  ispec.read(reinterpret_cast<char*>(&friction_perp_), sizeof(friction_perp_));
  ispec.read(reinterpret_cast<char*>(&bond_length_), sizeof(bond_length_));
  ispec.read(reinterpret_cast<char*>(&n_bonds_), sizeof(n_bonds_));
  bonds_.resize(n_bonds_, bonds_[0]);
  // Get initial site position
  for (int i=0; i<3; ++i)
    ispec.read(reinterpret_cast<char*>(&r0[i]), sizeof(double));
  for (int i=0; i<3; ++i)
    ispec.read(reinterpret_cast<char*>(&rf[i]), sizeof(double));
  // Initialize bonds from relative orientations
  for (int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    sites_[i_bond].SetPosition(r0);
    for (int i=0; i<3; ++i)
      ispec.read(reinterpret_cast<char*>(&u_bond[i]), sizeof(double));
    for (int i=0; i<n_dim_; ++i) {
      // Set bond position
      rf[i] = r0[i] + 0.5 * bond_length_ * u_bond[i];
      // Set next site position
      r0[i] += bond_length_ * u_bond[i];
    }
    sites_[i_bond].SetOrientation(u_bond);
    bonds_[i_bond].SetPosition(rf);
    bonds_[i_bond].SetOrientation(u_bond);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].SetLength(bond_length_);
    bonds_[i_bond].UpdatePeriodic();
  }
  sites_[n_bonds_].SetOrientation(sites_[n_bonds_-1].GetOrientation());
  double pos[3] = {0,0,0};
  //for (int i=0; i < params_->n_dim; ++i) {
    //pos[i] = bonds_[0].GetPosition()[i] - 0.5*bonds_[0].GetOrientation()[i] * bond_length_;
  //}
  //sites_[0].SetPosition(pos);
  for (int i=0; i < params_->n_dim; ++i) {
    pos[i] = bonds_[n_bonds_-1].GetPosition()[i] + 0.5*bonds_[n_bonds_-1].GetOrientation()[i] * bond_length_;
  }
  sites_[n_bonds_].SetPosition(pos);
  CalculateAngles();
}

void Filament::WritePosit(std::fstream &oposit) {
  double avg_pos[3], avg_u[3];
  GetAvgPosition(avg_pos);
  GetAvgOrientation(avg_u);
  std::copy(avg_pos,avg_pos+3,position_);
  UpdatePeriodic();
  for (auto& pos : position_)
    oposit.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for (auto& spos : scaled_position_)
    oposit.write(reinterpret_cast<char*>(&spos), sizeof(spos));
  for (auto& u : avg_u) 
    oposit.write(reinterpret_cast<char*>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char*>(&length_), sizeof(length_));
}

void Filament::ReadPosit(std::fstream &iposit) {
  if (iposit.eof()) return;
  double avg_pos[3], avg_u[3], s_pos[3];
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&avg_pos[i]), sizeof(double));
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&s_pos[i]), sizeof(double));
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&avg_u[i]), sizeof(double));
  iposit.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char*>(&length_), sizeof(length_));
  // Initialize first bond position
  for (int i=0; i<n_dim_; ++i) 
    avg_pos[i] = avg_pos[i] - 0.5*(length_ - bond_length_)*avg_u[i];
  for (int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    bonds_[i_bond].SetPosition(avg_pos);
    bonds_[i_bond].SetOrientation(avg_u);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].UpdatePeriodic();
    // Set next bond position
    for (int i=0; i<n_dim_; ++i) 
      avg_pos[i] += bond_length_*avg_u[i];
  }
  SetPosition(avg_pos);
  SetOrientation(avg_u);
  UpdatePeriodic();
}

void Filament::WriteCheckpoint(std::fstream &ocheck) {
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char*>(rng_state), rng_size);
  WriteSpec(ocheck);
}

void Filament::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char*>(rng_state), rng_size);
  ReadSpec(icheck);
  n_sites_ = n_bonds_+1;
  double r[3];
  for (int i=0; i<3; ++i)
    r[i] = bonds_[0].GetPosition()[i] - 0.5*bond_length_*bonds_[0].GetOrientation()[i];
  sites_[0].SetPosition(r);
  sites_[0].SetOrientation(bonds_[0].GetOrientation());
  for (int i_bond=1; i_bond<n_bonds_; ++i_bond) {
    for (int i=0; i<3; ++i)
      r[i] += bond_length_*sites_[i_bond-1].GetOrientation()[i];
    sites_[i_bond].SetPosition(r);
    sites_[i_bond].SetOrientation(bonds_[i_bond].GetOrientation());
  }
  for (int i=0; i<3; ++i)
    r[i] += bond_length_*sites_[n_sites_-2].GetOrientation()[i];
  sites_[n_sites_-1].SetPosition(r); 
  sites_[n_sites_-1].SetOrientation(bonds_[n_bonds_-1].GetOrientation());
  UpdateBondPositions();
  //Reallocate control structures
  tensions_.resize(n_sites_-1); //max_sites -1
  g_mat_lower_.resize(n_sites_-2); //max_sites-2
  g_mat_upper_.resize(n_sites_-2); //max_sites-2
  g_mat_diag_.resize(n_sites_-1); //max_sites-1
  det_t_mat_.resize(n_sites_+1); //max_sites+1
  det_b_mat_.resize(n_sites_+1); //max_sites+1
  g_mat_inverse_.resize(n_sites_-2); //max_sites-2
  k_eff_.resize(n_sites_-2); //max_sites-2
  h_mat_diag_.resize(n_sites_-1); //max_sites-1
  h_mat_upper_.resize(n_sites_-2); //max_sites-2
  h_mat_lower_.resize(n_sites_-2); //max_sites-2
  gamma_inverse_.resize(n_sites_*n_dim_*n_dim_); //max_sites*ndim*ndim
  cos_thetas_.resize(n_sites_-2); //max_sites-2
}

