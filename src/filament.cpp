#include <iostream>
#include <iomanip>
#include "filament.h"

void Filament::SetParameters(system_parameters *params) {
  length_ = params->rod_length;
  persistence_length_ = params->persistence_length;
  diameter_ = params->rod_diameter;
  // TODO JMM: add subdivisions of bonds for interactions, 
  //           should depend on cell length
  max_length_ = params->max_rod_length;
  min_length_ = params->min_rod_length;
  max_child_length_ = params->max_child_length;
  dynamic_instability_flag_ = params->dynamic_instability_flag;
  force_induced_catastrophe_flag_ = params->force_induced_catastrophe_flag;
  p_g2s_ = params->f_grow_to_shrink*delta_;
  p_g2p_ = params->f_grow_to_pause*delta_;
  p_s2p_ = params->f_shrink_to_pause*delta_;
  p_s2g_ = params->f_shrink_to_grow*delta_;
  p_p2s_ = params->f_pause_to_shrink*delta_;
  p_p2g_ = params->f_pause_to_grow*delta_;
  v_depoly_ = params->v_depoly;
  v_poly_ = params->v_poly;
  driving_factor_ = params->driving_factor;
  gamma_ratio_ = params->gamma_ratio;
  metric_forces_ = params->metric_forces;
  theta_validation_flag_ = params->theta_validation_flag;
  diffusion_validation_flag_ = params->diffusion_validation_flag;
}

void Filament::InitElements(system_parameters *params, space_struct *space) {
  n_bonds_ = (int) ceil(length_/max_child_length_);
  if (n_bonds_ < 2) 
    n_bonds_++;
  n_sites_ = n_bonds_+1;
  child_length_ = length_/n_bonds_;
  // If validating conformation of filament, create
  // the usual test filament
  if (theta_validation_flag_) {
    length_ = 8.0;
    max_child_length_ = 1.0;
    child_length_  = 1.0;
    min_length_ = 1.0;
    n_bonds_ = 8;
    dynamic_instability_flag_ = 0;
    force_induced_catastrophe_flag_ = 0;
  }
  if (length_/n_bonds_ < min_length_) {
    error_exit("ERROR: min_length_ of flexible filament segments too large for filament length.\n");
  }
  // Initialize sites
  for (int i=0; i<n_sites_; ++i) {
    Site s(params, space, gsl_rng_get(rng_.r), GetSID());
    s.SetCID(GetCID());
    elements_.push_back(s);
  }
  // Initialize bonds
  for (int i=0; i<n_bonds_; ++i) {
    Bond b(params, space, gsl_rng_get(rng_.r), GetSID());
    b.InitRID();
    b.SetCID(b.GetRID());
    //b.SetCID(GetCID());
    //b.SetRID(GetCID());
    v_elements_.push_back(b);
  }
  // Hack CIDs and RIDs to Set up interaction zones
  for (int i=1; i<n_bonds_-1; ++i) {
    if (i%2==0) 
      v_elements_[i].SetRID(v_elements_[i-1].GetRID());
    else
      v_elements_[i].SetCID(v_elements_[i-1].GetCID());
  }
  v_elements_[n_bonds_-1].SetCID(v_elements_[n_bonds_-2].GetCID());
  v_elements_[n_bonds_-1].SetRID(v_elements_[n_bonds_-2].GetRID());

  
  //Allocate control structures
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

void Filament::Init() {
  if (diffusion_validation_flag_) {
    DiffusionValidationInit();
    return;
  }
  InsertRandom();
  generate_random_unit_vector(n_dim_, orientation_, rng_.r);
  for (auto site=elements_.begin(); site!=elements_.end(); ++site) {
    site->SetDiameter(diameter_);
    site->SetLength(child_length_);
    site->SetPosition(position_);
    site->SetOrientation(orientation_);
    for (int i=0; i<n_dim_; ++i)
      position_[i] = position_[i] + orientation_[i] * child_length_;
    GenerateProbableOrientation();
  }
  UpdatePrevPositions();
  CalculateAngles();
  UpdateBondPositions();
  SetDiffusion();
  poly_state_ = GROW;
}

void Filament::DiffusionValidationInit() {
  for (int i=0; i<3; ++i) {
    position_[i] = 0.0;
    orientation_[i] = 0.0;
  }
  position_[n_dim_-1] = -0.5*length_;
  orientation_[n_dim_-1] = 1.0;
  for (auto site=elements_.begin(); site!=elements_.end(); ++site) {
    site->SetDiameter(diameter_);
    site->SetLength(child_length_);
    site->SetPosition(position_);
    site->SetOrientation(orientation_);
    for (int i=0; i<n_dim_; ++i)
      position_[i] = position_[i] + orientation_[i] * child_length_;
  }
  UpdatePrevPositions();
  CalculateAngles();
  UpdateBondPositions();
  SetDiffusion();
  poly_state_ = GROW;
}

void Filament::SetDiffusion() {
  double logLD = log(length_/diameter_);
  //double gamma_0 = 4.0/3.0*eps*((1+0.64*eps)/(1-1.15*eps) + 1.659 * SQR(eps));
  //gamma_perp_ = child_length_ * gamma_0;
  gamma_perp_ = 4.0*length_/(3.0*n_sites_*logLD);
  gamma_par_ = gamma_perp_ / gamma_ratio_;
  rand_sigma_perp_ = sqrt(24.0*gamma_perp_/delta_);
  rand_sigma_par_ = sqrt(24.0*gamma_par_/delta_);
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
    theta = acos( log( exp(-persistence_length_/child_length_) + 
          2.0*gsl_rng_uniform_pos(rng_.r)*sinh(persistence_length_/child_length_) ) 
        / (persistence_length_/child_length_) );
  }
  else {
    theta = acos( (log( 2.0*gsl_rng_uniform_pos(rng_.r)) - 
          log(2.0) + persistence_length_/child_length_)
        /(persistence_length_/child_length_) );
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
  rotate_orientation_vector(n_dim_, new_orientation, orientation_);
  std::copy(new_orientation,new_orientation+3,orientation_);
}

void Filament::UpdatePosition(bool midstep) {
  ApplyForcesTorques();
  Integrate(midstep);
  UpdateAvgPosition();
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
void Filament::Integrate(bool midstep) {
  CalculateAngles();
  CalculateTangents();
  if (midstep) {
    ConstructUnprojectedRandomForces();
    GeometricallyProjectRandomForces();
    UpdatePrevPositions();
  }
  AddRandomForces();
  CalculateBendingForces();
  CalculateTensions();
  UpdateSitePositions(midstep);
  UpdateBondPositions();
}

void Filament::CalculateAngles() {
  double cos_angle;
  for (int i_site=0; i_site<n_sites_-2; ++i_site) {
    double const * const u1 = elements_[i_site].GetOrientation();
    double const * const u2 = elements_[i_site+1].GetOrientation();
    cos_angle = dot_product(n_dim_, u1, u2);
    cos_thetas_[i_site] = cos_angle;
  }
}

void Filament::CalculateTangents() {
  elements_[0].SetTangent(elements_[0].GetOrientation());
  elements_[n_sites_-1].SetTangent(elements_[n_sites_-2].GetOrientation());
  double u_tan_mag;
  double u_tan[3];
  for (int i_site=1;i_site < n_sites_-1; ++i_site) {
    double const * const u1 = elements_[i_site-1].GetOrientation();
    double const * const u2 = elements_[i_site].GetOrientation();
    u_tan_mag = 0.0;
    for (int i=0; i<n_dim_; ++i) {
      u_tan[i] = u2[i]+u1[i];
      u_tan_mag += SQR(u_tan[i]);
    }
    u_tan_mag = sqrt(u_tan_mag);
    for (int i=0; i<n_dim_; ++i)
      u_tan[i]/=u_tan_mag;
    elements_[i_site].SetTangent(u_tan);
  }
}

void Filament::ConstructUnprojectedRandomForces() {
  // Create unprojected forces, see J. Chem. Phys. 122, 084903 (2005), eqn. 40.
  // xi is the random force vector with elements that are uncorrelated and randomly
  // distributed uniformly between -0.5 and 0.5, xi_term is the outer product of the
  // tangent vector u_tan_i u_tan_i acting on the vector xi
  double xi[3], xi_term[3], f_rand[3];
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    double const * const utan = elements_[i_site].GetTangent();
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
    elements_[i_site].SetRandomForce(f_rand);
  }
}

void Filament::GeometricallyProjectRandomForces() {
  double f_rand_temp[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    // Use the tensions vector to calculate the hard components of the random forces
    // These are not the same as the tensions, they will be calculated later
    double const * const f_rand1 = elements_[i_site].GetRandomForce();
    double const * const f_rand2 = elements_[i_site+1].GetRandomForce();
    double const * const u_site = elements_[i_site].GetOrientation();
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
    f_proj[i] = elements_[0].GetRandomForce()[i]
      + tensions_[0] * elements_[0].GetOrientation()[i];
  }
  elements_[0].SetRandomForce(f_proj);
  for (int i=0; i<n_dim_; ++i) {
    f_proj[i] = elements_[n_sites_-1].GetRandomForce()[i] 
    - tensions_[n_sites_-2] * elements_[n_sites_-2].GetOrientation()[i];
  }
  elements_[n_sites_-1].SetRandomForce(f_proj);
  // Then the rest
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    double const * const u1 = elements_[i_site-1].GetOrientation();
    double const * const u2 = elements_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      f_proj[i] = elements_[i_site].GetRandomForce()[i] + tensions_[i_site]*u2[i] - tensions_[i_site-1] * u1[i];
    }
    elements_[i_site].SetRandomForce(f_proj);
  }
}

void Filament::AddRandomForces() {
  for (auto site=elements_.begin(); site!=elements_.end(); ++site)
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
    k_eff_[i] = (persistence_length_ * (theta_validation_flag_ ? child_length_ : 1) + child_length_ * g_mat_inverse_[i])/SQR(child_length_);
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
        double const * const u1 = elements_[k_site-2].GetOrientation();
        double const * const u2 = elements_[k_site-1].GetOrientation();
        f_site[0] += k_eff_[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] );
        f_site[1] += k_eff_[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        double const * const u1 = elements_[k_site-1].GetOrientation();
        double const * const u2 = elements_[k_site].GetOrientation();
        f_site[0] += k_eff_[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1]
                        -((1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1]) );
        f_site[1] += k_eff_[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0]
                        -((1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0]) );
      }
      if (k_site<n_sites_-2) {
        double const * const u1 = elements_[k_site].GetOrientation();
        double const * const u2 = elements_[k_site+1].GetOrientation();
        f_site[0] -= k_eff_[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] );
        f_site[1] -= k_eff_[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0] );
      }
      elements_[k_site].AddForce(f_site);
    }
  }
  else if (n_dim_ == 3) {
    for (int k_site=0; k_site<n_sites_; ++k_site) {
      std::fill(f_site,f_site+3,0.0);
      if (k_site>1) {
        double const * const u1 = elements_[k_site-2].GetOrientation();
        double const * const u2 = elements_[k_site-1].GetOrientation();
        f_site[0] += k_eff_[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] );
        f_site[1] += k_eff_[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] );
        f_site[2] += k_eff_[k_site-2] * ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[2]*u2[1]*u1[1] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        double const * const u1 = elements_[k_site-1].GetOrientation();
        double const * const u2 = elements_[k_site].GetOrientation();
        f_site[0] += k_eff_[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] 
                        - ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] ) );
        f_site[1] += k_eff_[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2]
                        - ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] ) );
        f_site[2] += k_eff_[k_site-1] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1]
                        - ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[1]*u2[2]*u1[1] ) );
      }
      if(k_site<n_sites_-2) {
        double const * const u1 = elements_[k_site].GetOrientation();
        double const * const u2 = elements_[k_site+1].GetOrientation();
        f_site[0] -= k_eff_[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] );
        f_site[1] -= k_eff_[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2] );
        f_site[2] -= k_eff_[k_site] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1] );
      }
      elements_[k_site].AddForce(f_site);
    }
  }
}

void Filament::CalculateTensions() {
  // Calculate friction_inverse matrix
  int site_index = 0;
  int next_site = n_dim_*n_dim_;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    int gamma_index = 0;
    double const * const utan = elements_[i_site].GetTangent();
    for (int i=0; i<n_dim_; ++i) {
      for (int j=0; j<n_dim_; ++j) {
        gamma_inverse_[site_index+gamma_index] = 1.0/gamma_par_ * (utan[i]*utan[j])
                                     + 1.0/gamma_perp_ * ((i==j ? 1 : 0) - utan[i]*utan[j]);
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
    double const * const f1 = elements_[i_site].GetForce();
    double const * const f2 = elements_[i_site+1].GetForce();
    double const * const u2 = elements_[i_site].GetOrientation();
    double const * const utan1 = elements_[i_site].GetTangent();
    double const * const utan2 = elements_[i_site+1].GetTangent();
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
    h_mat_diag_[i_site] = 2.0/gamma_perp_ + (1.0/gamma_par_ - 1.0/gamma_perp_) *
      (SQR(utan1_dot_u2) + SQR(utan2_dot_u2));
    if (i_site>0) {
      double const * const u1 = elements_[i_site-1].GetOrientation();
      h_mat_upper_[i_site-1] = -1.0/gamma_perp_ * dot_product(n_dim_, u2, u1)
        - (1.0/gamma_par_ - 1.0/gamma_perp_) * (dot_product(n_dim_, utan1, u1) *
            dot_product(n_dim_, utan1, u2));
      h_mat_lower_[i_site-1] = h_mat_upper_[i_site-1];
    }
    site_index += next_site;
  }
  tridiagonal_solver(&h_mat_lower_, &h_mat_diag_, &h_mat_upper_, &tensions_, n_sites_-1);
}

void Filament::UpdateSitePositions(bool midstep) {
  double delta = (midstep ? 0.5*delta_ : delta_);
  double f_site[3];
  // First get total forces
  // Handle end sites first
  for (int i=0; i<n_dim_; ++i)
    f_site[i] = tensions_[0] * elements_[0].GetOrientation()[i];
  elements_[0].AddForce(f_site);
  for (int i=0; i<n_dim_; ++i)
    f_site[i] = -tensions_[n_sites_-2] * elements_[n_sites_-2].GetOrientation()[i];
  elements_[n_sites_-1].AddForce(f_site);
  // and then the rest
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    double const * const u_site1 = elements_[i_site-1].GetOrientation();
    double const * const u_site2 = elements_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      f_site[i] = tensions_[i_site] * u_site2[i] - tensions_[i_site-1] * u_site1[i];
    }
    elements_[i_site].AddForce(f_site);
  }
  // Now update positions
  double f_term[3], r_new[3];
  int site_index = 0;
  int next_site = n_dim_*n_dim_;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    double const * const f_site1 = elements_[i_site].GetForce();
    double const * const r_site1 = elements_[i_site].GetPosition();
    double const * const r_prev = elements_[i_site].GetPrevPosition();
    for (int i=0; i<n_dim_; ++i) {
      f_term[i] = gamma_inverse_[site_index+n_dim_*i] * f_site1[0] + gamma_inverse_[site_index+n_dim_*i+1] * f_site1[1];
      if (n_dim_ == 3)
        f_term[i] += gamma_inverse_[site_index+n_dim_*i+2] * f_site1[2];
      r_new[i] = r_prev[i] + delta * f_term[i];
    }
    elements_[i_site].SetPosition(r_new);
    site_index += next_site;
  }
  // Next, update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    double const * const r_site1 = elements_[i_site].GetPosition();
    double const * const r_site2 = elements_[i_site+1].GetPosition();
    u_mag = 0.0;
    for (int i=0; i<n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i=0; i<n_dim_; ++i)
      r_diff[i]/=u_mag;
    elements_[i_site].SetOrientation(r_diff);
  }
  elements_[n_sites_-1].SetOrientation(elements_[n_sites_-2].GetOrientation());
  // Finally, normalize site positions, making sure the sites are still rod-length apart
  for (int i_site=1; i_site<n_sites_; ++i_site) {
    double const * const r_site1 = elements_[i_site-1].GetPosition();
    double const * const u_site1 = elements_[i_site-1].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      r_diff[i] = r_site1[i] + child_length_ * u_site1[i];
    elements_[i_site].SetPosition(r_diff);
  }
}

void Filament::UpdateAvgPosition() {
  std::fill(position_, position_+3, 0.0);
  std::fill(orientation_, orientation_+3, 0.0);
  for (auto site_it : elements_) {
    double const * const site_pos = site_it.GetPosition();
    double const * const site_u = site_it.GetOrientation();
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

void Filament::UpdatePrevPositions() {
  for (auto site=elements_.begin(); site!=elements_.end(); ++site)
    site->SetPrevPosition(site->GetPosition());
}

double const * const Filament::GetDrTot() {return nullptr;}

void Filament::ApplyForcesTorques() {
  double pure_torque[3] = {0,0,0};
  double site_force[3] = {0,0,0};
  double linv=1.0/child_length_;
  for (int i=0; i<n_bonds_; ++i) {
    double const * const f = v_elements_[i].GetForce();
    double const * const t = v_elements_[i].GetTorque();
    double const * const u = elements_[i].GetOrientation();
    AddPotential(v_elements_[i].GetPotentialEnergy());
    // Convert torques into forces at bond ends
    // u x t / bond_length = pure torque force at tail of bond
    cross_product(u, t, pure_torque, 3);
    for (int i=0; i<n_dim_; ++i) {
      pure_torque[i]*=linv;
      site_force[i] = 0.5*f[i];
    }
    if (debug_trace) {
      //printf("potential %d: %2.8f\n",v_elements_[i].GetPotentialEnergy());
      //printf("pure torque %d: {%2.8f %2.8f %2.8f}\n",i,pure_torque[0],pure_torque[1],pure_torque[2]);
      //printf("site force %d: {%2.8f %2.8f %2.8f}\n",i,site_force[0],site_force[1],site_force[2]);
    }
    // Add translational forces and pure torque forces at bond ends
    elements_[i].AddForce(site_force);
    elements_[i].AddForce(pure_torque);
    for (int j=0; j<n_dim_; ++j)
      pure_torque[j] *= -1;
    elements_[i+1].AddForce(site_force);
    elements_[i+1].AddForce(pure_torque);
    // Add driving 
    double f_dr[3];
    for (int j=0; j<n_dim_; ++j)
      f_dr[j] = u[j]*driving_factor_;
    elements_[i].AddForce(f_dr);  

  }
}

void Filament::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetColor(color_, draw_type_);
    bond->Draw(graph_array);
  }
}

void Filament::UpdateBondPositions() {
  double pos[3];
  int i_site = 0;
  if (elements_.size() != v_elements_.size() + 1)
    error_exit("ERROR: n_bonds != n_sites - 1 in flexible filament!\n");
  for (auto bond=v_elements_.begin(); bond!=v_elements_.end(); ++bond) {
    double const * const r = elements_[i_site].GetPosition();
    double const * const u = elements_[i_site].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      pos[i] = r[i] + 0.5 * child_length_ * u[i];
    bond->SetPosition(pos);
    bond->SetOrientation(u);
    bond->SetLength(child_length_);
    bond->UpdatePeriodic();
    bond->SetRigidPosition(bond->GetPosition());
    bond->SetRigidScaledPosition(bond->GetScaledPosition());
    bond->SetRigidOrientation(u);
    bond->SetRigidLength(child_length_);
    bond->SetRigidDiameter(diameter_);
    i_site++;
  }
}

// Scale bond and site positions from new unit cell
void Filament::ScalePosition() {
  // scale first bond position using new unit cell
  v_elements_[0].ScalePosition();
  // then reposition sites based on first bond position
  // handle first site
  double r[3];
  double const * const bond_r = v_elements_[0].GetPosition();
  double const * const bond_u = v_elements_[0].GetOrientation();
  for (int i=0; i<n_dim_; ++i)
    r[i] = bond_r[i] - 0.5*child_length_*bond_u[i];
  elements_[0].SetPosition(r);
  // then handle remaining sites
  for (int i_site=1; i_site<n_sites_; ++i_site) {
    double const * const prev_r = elements_[i_site-1].GetPosition();
    double const * const prev_u = elements_[i_site-1].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      r[i] = prev_r[i] + child_length_*prev_u[i];
    elements_[i_site].SetPosition(r);
  }
  // update remaining bond positions
  UpdateBondPositions();
}

// Species specifics, including insertion routines
void FilamentSpecies::Configurator() {
  char *filename = params_->config_file;

  std::cout << "Filament species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";

  // Check insertion type
  std::string insertion_type;
  insertion_type = node["filament"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["filament"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap:        " << (can_overlap ? "true" : "false") << std::endl;

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 1; // default to orientation
  if (node["filament"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["filament"]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }
  if (node["filament"]["properties"]["draw_type"]) {
    std::string draw_type_s = node["filament"]["properties"]["draw_type"].as<std::string>();
    std::cout << "   draw_type: " << draw_type_s << std::endl;
    if (draw_type_s.compare("flat") == 0) {
      draw_type = 0;
    } else if (draw_type_s.compare("orientation") == 0) {
      draw_type = 1;
    }
  }

  if (insertion_type.compare("random") == 0) {
    int nfilaments          = node["filament"]["fil"]["num"].as<int>();
    double flength          = node["filament"]["fil"]["length"].as<double>();
    double plength          = node["filament"]["fil"]["persistence_length"].as<double>();
    double max_length       = node["filament"]["fil"]["max_length"].as<double>();
    double min_length       = node["filament"]["fil"]["min_length"].as<double>();
    double max_child_length = node["filament"]["fil"]["max_child_length"].as<double>();
    double diameter         = node["filament"]["fil"]["diameter"].as<double>();

    std::cout << std::setw(25) << std::left << "   n filaments:" << std::setw(10)
      << std::left << nfilaments << std::endl;
    std::cout << std::setw(25) << std::left << "   length:" << std::setw(10)
      << std::left << flength << std::endl;
    std::cout << std::setw(25) << std::left << "   persistence length:" << std::setw(10)
      << std::left << plength << std::endl;
    std::cout << std::setw(25) << std::left << "   max length:" << std::setw(10)
      << std::left << max_length << std::endl;
    std::cout << std::setw(25) << std::left << "   min length:" << std::setw(10)
      << std::left << min_length << std::endl;
    std::cout << std::setw(25) << std::left << "   max child length:" << std::setw(10)
      << std::left << max_child_length << std::endl;
    std::cout << std::setw(25) << std::left << "   diameter:" << std::setw(10)
      << std::left << diameter << std::endl;

    params_->n_filament = nfilaments;
    params_->rod_length = flength;
    params_->persistence_length = plength;
    params_->max_rod_length = max_length;
    params_->min_rod_length = min_length;
    params_->max_child_length = max_child_length; // Should depend on cell size
    params_->rod_diameter = diameter;

    n_members_ = 0;
    for (int i = 0; i < nfilaments; ++i) {
      Filament *member = new Filament(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init(); 

      if (can_overlap) {
        members_.push_back(member);
        n_members_++;
      } else {
        // Check for overlaps
        bool insert = true;
        auto new_simps = member->GetSimples();
        auto old_simps = GetSimples();
        for (auto ns : new_simps) {
          for (auto os : old_simps) {
            Interaction imd;
            MinimumDistance(ns, os, &imd, space_);
            if (imd.dr_mag2 < imd.buffer_mag2) {
              insert = false;
            }
          }
        }
        if (insert) {
          members_.push_back(member);
          n_members_++;
        }
        else {
          i--;
        }
      }
    }

  } else {
    // XXX Check for exact insertion of filaments
    std::cout << "Nope, not yet for filaments!\n";
    exit(1);
  }
  theta_validation_ = params_->theta_validation_flag ? true : false;
  n_dim_ = params_->n_dim;
  if (theta_validation_) {
    nbins_ = params_->n_bins;
    theta_distribution_ = new int**[n_members_];
    for (int i=0; i<n_members_; ++i) {
      theta_distribution_[i] = new int*[7]; // 7 angles to record for 8 bonds
      for (int j=0; j<7; ++j) {
        theta_distribution_[i][j] = new int[nbins_];
        std::fill(theta_distribution_[i][j],theta_distribution_[i][j]+nbins_,0.0);
      }
    }
  }
  midstep_ = true;
  ibin_ = 0;
  ivalidate_ = 0;
}

void FilamentSpecies::WriteThetaValidation(std::string run_name) {
  std::ostringstream file_name;
  file_name << run_name << "-thetas.log";
  std::ofstream thetas_file(file_name.str().c_str(), std::ios_base::out);
  thetas_file << "cos_theta" << " ";
  int n_theta = 1;
  for (int j_member=0; j_member<n_members_; ++j_member) {
    for (int k_angle=0; k_angle<7; ++k_angle) {
      thetas_file << "fil_" << j_member+1 << "_theta_" << n_theta++ << n_theta << " ";
    }
    n_theta = 1;
  }
  thetas_file << "\n";
  for (int i_bin=0; i_bin<nbins_; ++i_bin) {
    double angle = 2.0*i_bin/ (double) nbins_ - 1.0;
    thetas_file << angle << " ";
    for (int j_member=0; j_member<n_members_; ++j_member) {
      for (int k_angle=0; k_angle<7; ++k_angle) {
        thetas_file << theta_distribution_[j_member][k_angle][i_bin] << " ";
      }
    }
    thetas_file << "\n";
  }
}

void FilamentSpecies::WriteOutputs(std::string run_name) {
  if (theta_validation_) {
    WriteThetaValidation(run_name);
  }
}

void FilamentSpecies::ValidateThetaDistributions() {
  int i = 0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    std::vector<double> const * const thetas = (*it)->GetThetas();
    for (int j=0; j<7; ++j) {
      int bin_number = (int) floor( (1 + (*thetas)[j]) * (nbins_/2) );
      // Check boundaries
      if (bin_number == nbins_)
        bin_number = nbins_-1;
      else if (bin_number == -1)
        bin_number = 0;
      // Check for nonsensical values
      else if (bin_number > nbins_ || bin_number < 0) error_exit("Something went wrong in ValidateThetaDistributions!\n");
      theta_distribution_[i][j][bin_number]++;
    }
    i++;
  }
}

void Filament::GetAvgOrientation(double * au) {
  double avg_u[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=elements_.begin(); it!=elements_.end(); ++it) {
    double const * const u = it->GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      avg_u[i] += u[i];
    size++;
  }
  if (size == 0)
    error_exit("ERROR! Something went wrong in GetAvgOrientation!\n");
  for (int i=0; i<n_dim_; ++i)
    avg_u[i]/=size;
  std::copy(avg_u, avg_u+3, au);
}

void Filament::GetAvgPosition(double * ap) {
  double avg_p[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=elements_.begin(); it!=elements_.end(); ++it) {
    double const * const p = it->GetPosition();
    for (int i=0; i<n_dim_; ++i)
      avg_p[i] += p[i];
    size++;
  }
  if (size == 0)
    error_exit("ERROR! Something went wrong in GetAvgPosition!\n");
  for (int i=0; i<n_dim_; ++i)
    avg_p[i]/=size;
  std::copy(avg_p, avg_p+3, ap);
}

void Filament::DumpAll() {
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

/* The posit output for one filament is:
    diameter
    length
    bond_length
    n_bonds,
    position of first site
    position of last site
    all bond orientations
    */
void Filament::WritePosit(std::fstream &op){
  op.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  op.write(reinterpret_cast<char*>(&length_), sizeof(length_));
  op.write(reinterpret_cast<char*>(&child_length_), sizeof(child_length_));
  op.write(reinterpret_cast<char*>(&n_bonds_), sizeof(n_bonds_));
  double temp[3];
  double const * const r0 = elements_[0].GetPosition();
  std::copy(r0, r0+3, temp);
  for (auto& pos : temp)
    op.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  double const * const rf = elements_[n_bonds_].GetPosition();
  std::copy(rf, rf+3, temp);
  for (auto& pos : temp)
    op.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for (int i=0; i<n_bonds_; ++i) {
    double const * const orientation = v_elements_[i].GetOrientation();
    std::copy(orientation, orientation+3, temp);
    for (auto& u : temp) 
      op.write(reinterpret_cast<char*>(&u), sizeof(u));
  }
  return;
}

void Filament::ReadPosit(std::fstream &ip) {
  if (ip.eof()) return;
  double r0[3], rf[3], u_bond[3];
  Bond bond(params_, space_, params_->seed, sid_);
  ip.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ip.read(reinterpret_cast<char*>(&length_), sizeof(length_));
  ip.read(reinterpret_cast<char*>(&child_length_), sizeof(child_length_));
  ip.read(reinterpret_cast<char*>(&n_bonds_), sizeof(n_bonds_));
  v_elements_.resize(n_bonds_, bond);
  // Get initial site position
  for (auto& pos : r0)
    ip.read(reinterpret_cast<char*>(&pos), sizeof(pos));
  for (auto& pos : rf)
    ip.read(reinterpret_cast<char*>(&pos), sizeof(pos));
  // Initialize bonds from relative orientations
  for (int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    for (auto& u : u_bond)
      ip.read(reinterpret_cast<char*>(&u), sizeof(u));
    for (int i=0; i<n_dim_; ++i) {
      // Set bond position
      rf[i] = r0[i] + 0.5 * child_length_ * u_bond[i];
      // Set next site position
      r0[i] += child_length_ * u_bond[i];
    }
    v_elements_[i_bond].SetPosition(rf);
    v_elements_[i_bond].SetOrientation(u_bond);
    v_elements_[i_bond].SetDiameter(diameter_);
    v_elements_[i_bond].SetLength(child_length_);
    v_elements_[i_bond].UpdatePeriodic();
  }
}

