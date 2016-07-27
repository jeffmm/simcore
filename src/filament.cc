#include <iostream>
#include "filament.h"

// Constructor for array reference construction
Filament::Filament() {}

// Initialization constructor
Filament::Filament(system_parameters *params, SpaceProperties *space, control_structure *control, 
    attach_state_t attach_state, long seed) {
 
  space_ = space;
  params_ = params;
  ctrl_ = control;
  length_ = 0;
  diameter_ = params_->filament_diameter;
  attach_state_ = attach_state;
  poly_state_ = GROW;
  int max_sites = control->max_sites;
  sites.reserve(max_sites);
  cos_thetas_ = new double[max_sites-2];
  a_rand_ = sqrt(24.0/params_->delta);
  midstep_ = true;
  rng_.init(seed);
  first_iteration_ = true;
}

void Filament::InitRandomFree() {
  int n_dim = params_->n_dim;
  int i_segment;
  bool in_bounds = false;
  double *trial_orientation = new double[n_dim];
  double *trial_coord;
  do {
    length_ = params_->min_length + (params_->max_length-params_->min_length) * gsl_rng_uniform_pos(rng_.r);
    n_segments_ = (int) ceil(length_/params_->max_segment_length);
    n_sites_ = n_segments_ + 1;
    segment_length_ = length_/(n_segments_);
    SetDiffusion();
    i_segment = 0;
    sites.clear();
    trial_coord = space_->RandomCoordinate(0.5*diameter_);
    generate_random_unit_vector(n_dim, trial_orientation, rng_.r);
    Site init_site(n_dim);
    init_site.SetPosition(trial_coord);
    init_site.SetScaledPosition(trial_coord);
    init_site.SetOrientation(trial_orientation);
    init_site.SetLength(segment_length_);
    init_site.SetDiameter(diameter_);
    sites.push_back(init_site);
    in_bounds = true;
    do {
      for (int i=0; i<n_dim; ++i) {
        trial_coord[i] = trial_coord[i] + trial_orientation[i] * segment_length_;
      }
      trial_orientation = GenerateProbableOrientation(trial_orientation);
      Site trial_site(n_dim);
      trial_site.SetPosition(trial_coord);
      trial_site.SetScaledPosition(trial_coord);
      trial_site.SetOrientation(trial_orientation);
      trial_site.SetLength(segment_length_);
      trial_site.SetDiameter(diameter_);
      sites.push_back(trial_site);
      i_segment++;
      in_bounds = space_->CheckInBounds(trial_coord, 0.5*diameter_);
      // If using snowman boundary, it's not enough to make sure the beads are
      // in bounds to check if we have violations of the boundary condition,
      // due to the concavity of the snowman boundary, must check segments
      if (space_->GetType() == SNOWMAN && in_bounds) {
        int ind = (int) sites.size() - 1;
        in_bounds = space_->CheckSegmentInBounds(sites[ind].GetPosition(), sites[ind-1].GetPosition(), 
            0.5*diameter_);
      }
    } while (i_segment < n_segments_ && in_bounds);
  } while (!in_bounds);
  UpdatePrevPositions();
  CalculateAngles(false);
  ZeroForces();
}

void Filament::InitRandomAttached(double *sphere_position, double sphere_radius) {
  int n_dim = params_->n_dim;
  int i_segment;
  bool in_bounds = false;
  double *trial_orientation = new double[n_dim];
  double *trial_attach = new double[n_dim];
  double *trial_coord = new double[n_dim];
  do {
    length_ = params_->min_length + (params_->max_length-params_->min_length) * gsl_rng_uniform_pos(rng_.r);
    n_segments_ = (int) ceil(length_/params_->max_segment_length);
    n_sites_ = n_segments_ + 1;
    segment_length_ = length_/(n_segments_);
    SetDiffusion();
    i_segment = 0;
    sites.clear();
    generate_random_unit_vector(n_dim, trial_attach, rng_.r);
    trial_orientation = GenerateProbableOrientation(trial_attach);
    for (int i=0; i<n_dim; ++i) {
      anchor_.r_anchor[i] = sphere_position[i] + sphere_radius * trial_attach[i];
      anchor_.u_anchor[i] = trial_attach[i];
      trial_coord[i] = anchor_.r_anchor[i] + 1e-8 * trial_attach[i];
    }
    Site init_site(n_dim);
    init_site.SetPosition(trial_coord);
    init_site.SetScaledPosition(trial_coord);
    init_site.SetOrientation(trial_orientation);
    init_site.SetLength(segment_length_);
    init_site.SetDiameter(diameter_);
    sites.push_back(init_site);
    in_bounds = true;
    do {
      for (int i=0; i<n_dim; ++i) {
        trial_coord[i] = trial_coord[i] + trial_orientation[i] * segment_length_;
      }
      trial_orientation = GenerateProbableOrientation(trial_orientation);
      Site trial_site(n_dim);
      trial_site.SetPosition(trial_coord);
      trial_site.SetScaledPosition(trial_coord);
      trial_site.SetOrientation(trial_orientation);
      trial_site.SetLength(segment_length_);
      trial_site.SetDiameter(diameter_);
      sites.push_back(trial_site);
      i_segment++;
      in_bounds = space_->CheckInBounds(trial_coord, 0.5*diameter_);
      // If using snowman boundary, it's not enough to make sure the beads are
      // in bounds to check if we have violations of the boundary condition,
      // due to the concavity of the snowman boundary, must check segments
      if (space_->GetType() == SNOWMAN && in_bounds) {
        int ind = (int) sites.size() - 1;
        in_bounds = space_->CheckSegmentInBounds(sites[ind].GetPosition(), sites[ind-1].GetPosition(),
            0.5*diameter_);
      }
    } while (i_segment < n_segments_ && in_bounds);
  } while (!in_bounds);
  CalculateAngles(false);
  ZeroForces();
}

void Filament::SetDiffusion() {
  double eps = log(2.0*length_);
  double friction_0 = 4.0/3.0*eps*((1+0.64*eps)/(1-1.15*eps) + 1.659 * SQR(eps));
  friction_perp_ = segment_length_ * friction_0;
  friction_par_ = friction_perp_ / params_->friction_ratio;
  sqrt_friction_perp_ = sqrt(friction_perp_);
  sqrt_friction_par_ = sqrt(friction_par_);
}
  

double *Filament::GenerateProbableOrientation(double *vec) {

  // Generate random theta pulled from probability distribution P(th) = exp(k cos(th))
  // where k is the persistence_length of the filament
  // If k is too large, there is enormous imprecision in this calculation since sinh(k) is very large
  // so to fix this I introduce an approximate distribution that is valid for large k
  double theta;
  int n_dim = params_->n_dim;

  if (params_->persistence_length == 0) {
    theta = gsl_rng_uniform_pos(rng_.r) * M_PI;
  }
  else if (params_->persistence_length < 100) {
    theta = acos( log( exp(-params_->persistence_length/segment_length_) + 
          2.0*gsl_rng_uniform_pos(rng_.r)*sinh(params_->persistence_length/segment_length_) ) 
        / (params_->persistence_length/(params_->theta_validation_flag ? 1 : segment_length_)) );
  }
  else {
    theta = acos( (log( 2.0*gsl_rng_uniform_pos(rng_.r)) - 
          log(2.0) + params_->persistence_length/segment_length_)
        /(params_->persistence_length/(params_->theta_validation_flag ? 1 : segment_length_)) );
  }

  double *new_vec = new double[n_dim];
  if (n_dim==2) {
    theta = (gsl_rng_uniform_int(rng_.r,2)==0 ? -1 : 1) * theta;
    new_vec[0] = cos(theta);
    new_vec[1] = sin(theta);
  }
  else {
    double phi = gsl_rng_uniform_pos(rng_.r) * 2.0 * M_PI;
    new_vec[0] = sin(theta)*cos(phi);
    new_vec[1] = sin(theta)*sin(phi);
    new_vec[2] = cos(theta);
  }
  rotate_orientation_vector(n_dim, new_vec, vec);
  return new_vec;
}

void Filament::UpdatePosition() {
  CalculateAngles(false);
  CalculateTangents();
  if (midstep_) {
    if (params_->theta_validation_flag) {
      ValidateThetaDistribution();
    }
    if (params_->position_correlation_flag) {
      ValidatePositionStepDistribution();
    }
    GenerateRandomForces();
    ProjectRandomForces();
    UpdatePrevPositions();
  }
  AddRandomForces();
  CalculateBendingForces();
  CalculateTensions();
  UpdateSitePositions();
  midstep_ = !midstep_;
}

void Filament::CalculateTangents() {
  int n_dim = params_->n_dim;
  double *utan1 = sites[0].GetTangent();
  double *u1 = sites[0].GetOrientation();
  double *utan2 = sites[n_sites_-1].GetTangent();
  double *u2 = sites[n_sites_-2].GetOrientation();
  for (int i=0; i<n_dim; ++i) {
    utan1[i] = u1[i];
    utan2[i] = u2[i];
  }
  double u_mag;
  double u_add[3];
  double *utan;
  u2=u1;
  for (int i_site=1;i_site < n_sites_-1; ++i_site) {
    u1=u2;
    u2 = sites[i_site].GetOrientation();
    u_mag = 0.0;
    for (int i=0; i<n_dim; ++i) {
      u_add[i] = u2[i]+u1[i];
      u_mag += SQR(u_add[i]);
    }
    u_mag = sqrt(u_mag);
    utan = sites[i_site].GetTangent();
    for (int i=0; i<n_dim; ++i)
      utan[i] = u_add[i]/u_mag;
  }
}

void Filament::CalculateAngles(bool rescale) {
  double *u1, *u2;
  int n_dim = params_->n_dim;
  double cos_angle;
  bool sharp_angle = false;
  for (int i_site=0; i_site<n_sites_-2; ++i_site) {
    u1 = sites[i_site].GetOrientation();
    u2 = sites[i_site+1].GetOrientation();
    cos_angle = dot_product(n_dim, u1, u2);
    cos_thetas_[i_site] = cos_angle;
    if (cos_angle <= 0 && params_->dynamic_instability_flag)
      error_exit("ERROR: Acute angle between adjoining segments detected with dynamic instabiliy on.\n \
          Increase persistence length or turn off dynamic instability for floppy filaments\n");
    if (cos_angle < 0.7 && params_->dynamic_instability_flag)
      sharp_angle = true;
  }
  // If dynamic instability is on, make sure the angle between two adjoining segments is
  // less than 0.5*pi, or else the site rescaling will fail horribly.
  // This is only an issue for floppy filaments.
  if (rescale && sharp_angle && midstep_ && length_/(n_segments_+1) > params_->min_segment_length) {
    AddSite();
    CalculateAngles(true);
  }
}


void Filament::UpdatePrevPositions() {
  for (site_iterator i_site = sites.begin(); i_site != sites.end(); ++i_site) 
    i_site->SetPrevPosition();
}

void Filament::AddRandomForces() {
  for (site_iterator i_site = sites.begin(); i_site != sites.end(); ++i_site) 
    i_site->AddRandomForce();
}

void Filament::ValidateThetaDistribution() {
  int n_bins = params_->n_bins;
  int bin_number;
  int **cos_dist = ctrl_->cos_theta_dist;
  for (int i=0; i<n_segments_-1; ++i) {
    bin_number = (int) floor( (1 + cos_thetas_[i]) * (n_bins/2) );
    if (bin_number == n_bins) bin_number = n_bins-1;
    else if (bin_number == -1) bin_number = 0;
    else if (bin_number > n_bins && bin_number < 0) error_exit("Something went wrong in validate_theta_dist_sphero!\n");
    cos_dist[i][bin_number]++;
  }
}

void Filament::ValidatePositionStepDistribution() {
  if (first_iteration_) {
    first_iteration_ = false;
    return; // FIXME
  }
  int n_dim = params_->n_dim;
  double *dr = new double[n_dim];
  double *r_old, *r_new;
  int **pos_corr = ctrl_->position_step_correlation;
  int n_site = 0;
  for (site_iterator i_site = sites.begin(); i_site != sites.end(); ++i_site) {
    r_old = i_site->GetPrevPosition();
    r_new = i_site->GetPosition();
    for (int i=0; i<n_dim; ++i)
      dr[i] = r_new[i] - r_old[i];
    normalize_vector(dr, n_dim);
    double correlation = dot_product(n_dim, i_site->GetOrientation(), dr);
    int bin_number = (int) floor( (1 + correlation) * (params_->n_bins/2) );
    pos_corr[n_site][bin_number]++;
    n_site++;
  }
}

void Filament::GenerateRandomForces() {
  int n_dim = params_->n_dim;
  double xi[3];
  double xi_term[3];
  double *utan, *f_rand;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    utan = sites[i_site].GetTangent();
    f_rand = sites[i_site].GetRandomForce();
    for (int i=0; i<n_dim; ++i) 
      xi[i] = gsl_rng_uniform_pos(rng_.r) - 0.5;
    // Create unprojected forces, see J. Chem. Phys. 122, 084903 (2005), eqn. 40.
    if (n_dim == 2) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0];
    }
    else if (n_dim==3) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1]
                  + utan[0] * utan[2] * xi[2];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0]
                  + utan[1] * utan[2] * xi[2];
      xi_term[2] = SQR(utan[2]) * xi[2] + utan[0] * utan[2] * xi[0]
                  + utan[1] * utan[2] * xi[1];
    }
    for (int i=0; i<n_dim; ++i) {
      f_rand[i] = a_rand_ * (sqrt_friction_perp_ * xi[i] + 
          (sqrt_friction_par_ - sqrt_friction_perp_) * xi_term[i]);
    }
  }
}

void Filament::ProjectRandomForces() {

  int n_dim = params_->n_dim;
  double *p_vec = ctrl_->p_vec;
  double *g_mat_upper = ctrl_->g_mat_upper;
  double *g_mat_lower = ctrl_->g_mat_lower;
  double *g_mat_diag = ctrl_->g_mat_diag;

  double f_rand_temp[3];
  double *f_rand1, *f_rand2, *u_site;
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    // First get the p_vec elements
    f_rand1 = sites[i_site].GetRandomForce();
    f_rand2 = sites[i_site+1].GetRandomForce();
    u_site = sites[i_site].GetOrientation();
    for (int i=0; i<n_dim; ++i)
      f_rand_temp[i] = f_rand2[i] - f_rand1[i];
    p_vec[i_site] = dot_product(n_dim, f_rand_temp, u_site);
    // Then get the G arrays (for inertialess case where m=1, see
    // ref. 15 of above paper)
    g_mat_diag[i_site] = 2;
    if (i_site > 0) {
      g_mat_upper[i_site-1] = - cos_thetas_[i_site-1];
      g_mat_lower[i_site-1] = - cos_thetas_[i_site-1];
    }
  }

  // Now solve using tridiagonal solver
  tridiagonal_solver(g_mat_lower, g_mat_diag, g_mat_upper, p_vec, n_sites_-1);

  // Update to the projected brownian forces
  // First the end sites:
  
  f_rand1 = sites[0].GetRandomForce();
  f_rand2 = sites[n_sites_-1].GetRandomForce();
  u_site = sites[0].GetOrientation();
  double *u_site2 = sites[n_sites_-2].GetOrientation();
  for (int i=0; i<n_dim; ++i) {
    f_rand1[i] += p_vec[0] * u_site[i];
    f_rand2[i] -= p_vec[n_sites_-2] * u_site2[i];
  }
  // Then the rest
  u_site2=u_site;
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    f_rand1 = sites[i_site].GetRandomForce();
    u_site=u_site2;
    u_site2 = sites[i_site].GetOrientation();
    for (int i=0; i<n_dim; ++i) {
      f_rand1[i] +=  p_vec[i_site]*u_site2[i] - p_vec[i_site-1] * u_site[i];
    }
  }
}

void Filament::CalculateBendingForces() {

  int n_dim = params_->n_dim;
  double *det_t_mat = ctrl_->det_t_mat;
  double *det_b_mat = ctrl_->det_b_mat;
  double *g_mat_upper = cos_thetas_;
  double *g_mat_inverse = ctrl_->g_mat_inverse;
  double *k_eff = ctrl_->k_eff;

  if (params_->metric_forces) {
    det_t_mat[0] = 1;
    det_t_mat[1] = 2;
    det_b_mat[n_sites_] = 1;
    det_b_mat[n_sites_-1] = 2;
    for (int i=2; i<n_sites_; ++i) {
      det_t_mat[i] = 2 * det_t_mat[i-1] - SQR(- g_mat_upper[i-2]) * det_t_mat[i-2];
      det_b_mat[n_sites_-i] = 2 * det_b_mat[n_sites_-i+1] - SQR(- g_mat_upper[n_sites_-i-1]) * det_b_mat[n_sites_-i+2];
    }
    double det_g = det_t_mat[n_sites_-1];
    for(int i=0; i<n_sites_-2; ++i) {
      g_mat_inverse[i] = g_mat_upper[i] * det_t_mat[i] * det_b_mat[i+3] / det_g;
    }
  }
  else {
    for(int i=0; i<n_sites_-2; ++i) {
      g_mat_inverse[i] = 0;
    }
  }
  // Now calculate the effective rigidities
  for (int i=0; i<n_sites_-2; ++i) {
    k_eff[i] = (params_->persistence_length * ( params_->theta_validation_flag ? segment_length_ : 1) + segment_length_ * g_mat_inverse[i])/SQR(segment_length_);
  }
  // Using these, we can calculate the forces on each of the sites

  // These calculations were done by hand and are not particularly readable,
  // but are more efficient than doing it explicitly in the code for readability
  // If this ever needs fixed, you need to either check the indices very carefully
  // or redo the calculation by hand! 
  // See Pasquali and Morse, J. Chem. Phys. Vol 116, No 5 (2002)
  double *u1, *u2, *f_site;
  if (n_dim == 2) {
    for (int k_site=0; k_site<n_sites_; ++k_site) {
      f_site = sites[k_site].GetForce();
      if (k_site>1) {
        u1 = sites[k_site-2].GetOrientation();
        u2 = sites[k_site-1].GetOrientation();
        f_site[0] += k_eff[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] );
        f_site[1] += k_eff[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        u1 = sites[k_site-1].GetOrientation();
        u2 = sites[k_site].GetOrientation();
        f_site[0] += k_eff[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1]
                        -((1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1]) );
        f_site[1] += k_eff[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0]
                        -((1-SQR(u2[1]))*u1[1] - u2[0]*u2[1]*u1[0]) );
      }
      if (k_site<n_sites_-2) {
        u1 = sites[k_site].GetOrientation();
        u2 = sites[k_site+1].GetOrientation();
        f_site[0] -= k_eff[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] );
        f_site[1] -= k_eff[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[0]*u1[1]*u2[0] );
      }
    }
  }
  else if (n_dim == 3) {
    for (int k_site=0; k_site<n_sites_; ++k_site) {
      f_site = sites[k_site].GetForce();
      if (k_site>1) {
        u1 = sites[k_site-2].GetOrientation();
        u2 = sites[k_site-1].GetOrientation();
        f_site[0] += k_eff[k_site-2] * ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] );
        f_site[1] += k_eff[k_site-2] * ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] );
        f_site[2] += k_eff[k_site-2] * ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[2]*u2[1]*u1[1] );
      }
      if (k_site>0 && k_site<n_sites_-1) {
        u1 = sites[k_site-1].GetOrientation();
        u2 = sites[k_site].GetOrientation();
        f_site[0] += k_eff[k_site-1] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] 
                        - ( (1-SQR(u2[0]))*u1[0] - u2[0]*u2[1]*u1[1] - u2[0]*u2[2]*u1[2] ) );
        f_site[1] += k_eff[k_site-1] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2]
                        - ( (1-SQR(u2[1]))*u1[1] - u2[1]*u2[0]*u1[0] - u2[1]*u2[2]*u1[2] ) );
        f_site[2] += k_eff[k_site-1] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1]
                        - ( (1-SQR(u2[2]))*u1[2] - u2[2]*u2[0]*u1[0] - u2[1]*u2[2]*u1[1] ) );
      }
      if(k_site<n_sites_-2) {
        u1 = sites[k_site].GetOrientation();
        u2 = sites[k_site+1].GetOrientation();
        f_site[0] -= k_eff[k_site] * ( (1-SQR(u1[0]))*u2[0] - u1[0]*u1[1]*u2[1] - u1[0]*u1[2]*u2[2] );
        f_site[1] -= k_eff[k_site] * ( (1-SQR(u1[1]))*u2[1] - u1[1]*u1[0]*u2[0] - u1[1]*u1[2]*u2[2] );
        f_site[2] -= k_eff[k_site] * ( (1-SQR(u1[2]))*u2[2] - u1[2]*u1[0]*u2[0] - u1[2]*u1[1]*u2[1] );
      }
    }
  }
}

void Filament::CalculateTensions() {

  int n_dim = params_->n_dim;
  double **friction_inverse = ctrl_->friction_inverse;
  double *h_mat_upper = ctrl_->h_mat_upper;
  double *h_mat_lower = ctrl_->h_mat_lower;
  double *h_mat_diag = ctrl_->h_mat_diag;
  double *tensions = ctrl_->p_vec;

  // Calculate friction_inverse matrix
  double *utan;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    int ind=0;
    utan = sites[i_site].GetTangent();
    for (int i=0; i<n_dim; ++i) {
      for (int j=0; j<n_dim; ++j) {
        friction_inverse[i_site][ind] = 1.0/friction_par_ * (utan[i]*utan[j])
                                     + 1.0/friction_perp_ * ((i==j ? 1 : 0) - utan[i]*utan[j]);
        ind++;
      }
    }
  }
  // Populate the H matrix and Q vector (using tensions array)
  double temp_a, temp_b;
  double f_diff[3];
  double *f1, *f2, *u1, *u2, *utan1, *utan2;
  double utan1_dot_u2, utan2_dot_u2;
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    // f_diff is the term in par_entheses in equation 29 of J. Chem. Phys. 122, 084903 (2005)
    f1 = sites[i_site].GetForce();
    f2 = sites[i_site+1].GetForce();
    u2 = sites[i_site].GetOrientation();
    utan1 = sites[i_site].GetTangent();
    utan2 = sites[i_site+1].GetTangent();
    for (int i=0; i<n_dim; ++i) {
      temp_a = friction_inverse[i_site][n_dim*i] * f1[0] + friction_inverse[i_site][n_dim*i+1] * f1[1];
      if (n_dim == 3) temp_a += friction_inverse[i_site][n_dim*i+2] * f1[2];
      temp_b = friction_inverse[i_site+1][n_dim*i] * f2[0] + friction_inverse[i_site+1][n_dim*i+1] * f2[1];
      if (n_dim == 3) temp_b += friction_inverse[i_site+1][n_dim*i+2] * f2[2];
      f_diff[i] = temp_b - temp_a;
    }
    tensions[i_site] = dot_product(n_dim, u2, f_diff);
    utan1_dot_u2 = dot_product(n_dim, utan1, u2);
    utan2_dot_u2 = dot_product(n_dim, utan2, u2);
    h_mat_diag[i_site] = 2.0/friction_perp_ + (1.0/friction_par_ - 1.0/friction_perp_) *
      (SQR(utan1_dot_u2) + SQR(utan2_dot_u2));
    if (i_site>0) {
      u1 = sites[i_site-1].GetOrientation();
      h_mat_upper[i_site-1] = -1.0/friction_perp_ * dot_product(n_dim, u2, u1)
        - (1.0/friction_par_ - 1.0/friction_perp_) * (dot_product(n_dim, utan1, u1) *
            dot_product(n_dim, utan1, u2));
      h_mat_lower[i_site-1] = h_mat_upper[i_site-1];
    }
  }
  tridiagonal_solver(h_mat_lower, h_mat_diag, h_mat_upper, tensions, n_sites_-1);
}

void Filament::UpdateSitePositions() {

  int n_dim = params_->n_dim;
  double delta = (midstep_ ? 0.5*params_->delta : params_->delta);
  double *tensions = ctrl_->p_vec;
  double **friction_inverse = ctrl_->friction_inverse;

  // Get total forces
  // First get the end sites
  double *f_site1, *f_site2, *u_site1, *u_site2, *r_site1, *r_site2, *r_prev;
  f_site1 = sites[0].GetForce();
  f_site2 = sites[n_sites_-1].GetForce();
  u_site1 = sites[0].GetOrientation();
  u_site2 = sites[n_sites_-2].GetOrientation();
  for (int i=0; i<n_dim; ++i) {
    f_site1[i] += tensions[0] * u_site1[i];
    f_site2[i] -= tensions[n_sites_-2] * u_site2[i];
  }
  // and then the rest
  for (int i_site=1; i_site<n_sites_-1; ++i_site) {
    f_site1 = sites[i_site].GetForce();
    u_site1 = sites[i_site-1].GetOrientation();
    u_site2 = sites[i_site].GetOrientation();
    for (int i=0; i<n_dim; ++i) {
      f_site1[i] += tensions[i_site] * u_site2[i] - tensions[i_site-1] * u_site1[i];
    }
  }

  // Now update positions
  double f_term[3];
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    f_site1 = sites[i_site].GetForce();
    r_site1 = sites[i_site].GetPosition();
    r_prev = sites[i_site].GetPrevPosition();
    for (int i=0; i<n_dim; ++i) {
      f_term[i] = friction_inverse[i_site][n_dim*i] * f_site1[0] + friction_inverse[i_site][n_dim*i+1] * f_site1[1];
      if (n_dim == 3) f_term[i] += friction_inverse[i_site][n_dim*i+2] * f_site1[2];
      r_site1[i] = r_prev[i] + delta * f_term[i];
    }
  }

  // Update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    r_site1 = sites[i_site].GetPosition();
    r_site2 = sites[i_site+1].GetPosition();
    u_mag = 0.0;
    for (int i=0; i<n_dim; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    if (params_->error_analysis_flag) {
      ctrl_->pos_error[ctrl_->pos_error_ind] = ABS(segment_length_ - u_mag)/segment_length_;
      ctrl_->pos_error_ind++;
    }
    u_site1 = sites[i_site].GetOrientation();
    for (int i=0; i<n_dim; ++i)
      u_site1[i] = (r_diff[i])/u_mag;
  }
  sites[n_sites_-1].SetOrientation(sites[n_sites_-2].GetOrientation());
  // Normalize site positions... in other words, make sure the sites are still rod-length apart
  r_site2 = sites[0].GetPosition();
  for (int i_site=1; i_site<n_sites_; ++i_site) {
    r_site1 = r_site2;
    r_site2 = sites[i_site].GetPosition();
    u_site1 = sites[i_site-1].GetOrientation();
    for (int i=0; i<n_dim; ++i)
      r_site2[i] = r_site1[i] + segment_length_ * u_site1[i];
  }
}

void Filament::DynamicInstability() {
  UpdatePolyState();
  GrowFilament();
  SetDiffusion();
}

void Filament::GrowFilament() {
  // If the filament is paused, do nothing
  if (poly_state_ == PAUSE)
    return;
  double delta = params_->delta;
  double delta_length;
  if (poly_state_ == GROW) {
    delta_length = params_->v_poly * delta;
    length_ += delta_length;
  }
  else if (poly_state_ == SHRINK) {
    delta_length = params_->v_depoly * delta;
    length_ -= delta_length;
  }
  RescaleSegments();
  if (segment_length_ > params_->max_segment_length) 
    AddSite();
  else if (segment_length_ < params_->min_segment_length) 
    RemoveSite();
}

void Filament::RescaleSegments() {
  UpdatePrevPositions();
  int n_dim = params_->n_dim;
  double old_segment_length = segment_length_;
  segment_length_ = length_/n_segments_;
  double k, dl;
  if (poly_state_ == SHRINK) {
    double *r2=sites[1].GetPosition();
    double *r1=sites[0].GetPrevPosition();
    double *u1=sites[0].GetOrientation();
    double *r2_old;
    for (int i=0; i<n_dim; ++i) {
      r2[i] = r1[i] + u1[i] * segment_length_;
    }
    dl = old_segment_length - segment_length_;
    for (int i_site=2; i_site<n_sites_; ++i_site) {
      r2 = sites[i_site].GetPosition();
      r2_old = sites[i_site].GetPrevPosition();
      r1 = sites[i_site-1].GetPrevPosition();
      u1 = sites[i_site-1].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-2]) - SQR(dl) + SQR(segment_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site-2]*dl + sqrt(k);
      for (int i=0; i<n_dim; ++i)
        r2[i] = r1[i] + k * u1[i];
      dl = old_segment_length - k;
    }
  }
  else if (poly_state_ == GROW) {
    double *r, *r_old, *u;
    dl = old_segment_length;
    for (int i_site=1; i_site<n_sites_-1; ++i_site) {
      r = sites[i_site].GetPosition();
      r_old = sites[i_site].GetPrevPosition();
      u = sites[i_site].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-1]) - SQR(dl) + SQR(segment_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site-1]*dl + sqrt(k);
      for (int i=0; i<n_dim; ++i)
        r[i] = r_old[i] + k * u[i];
      dl = old_segment_length - k;
    }
    r = sites[n_sites_-1].GetPosition();
    u = sites[n_sites_-2].GetOrientation();
    r_old = sites[n_sites_-2].GetPosition();
    for (int i=0; i<n_dim; ++i)
      r[i] = r_old[i] + segment_length_ * u[i];
  }
  for (site_iterator i_site=sites.begin();
      i_site != sites.end();
      ++i_site)
    i_site->SetLength(segment_length_);
  UpdateSiteOrientations();
  CalculateAngles(false);
}

void Filament::UpdateSiteOrientations() {
  int n_dim = params_->n_dim;
  for (site_iterator i_site = sites.begin();
      i_site != sites.end()-1;
      ++i_site) {
    double u_site[3];
    double *r_site1 = i_site->GetPosition();
    double *r_site2 = (i_site+1)->GetPosition();
    for (int i=0; i<n_dim; ++i)
      u_site[i] = r_site2[i] - r_site1[i];
    normalize_vector(u_site, n_dim);
    i_site->SetOrientation(u_site);
  }
  // Have the orientation of the last site be parallel to the orientation
  // vector of the segment before it
  sites[n_sites_-1].SetOrientation(sites[n_sites_-2].GetOrientation());
}

void Filament::AddSite() {
  UpdatePrevPositions();
  int n_dim = params_->n_dim;
  Site new_site(n_dim);
  new_site.SetDiameter(diameter_);
  // Place the new site in the same position as the current last site
  sites.push_back(new_site);
  n_sites_++;
  n_segments_++;
  double old_segment_length = segment_length_;
  segment_length_ = length_/n_segments_;
  double error = 0;
  double dl = 1;
  int j=0;
  double *r2, *r1, *u1, *r2_old, k;
  do {
    segment_length_ = segment_length_ + error/n_segments_;
    r2=sites[1].GetPosition();
    r1=sites[0].GetPrevPosition();
    u1=sites[0].GetOrientation();
    r2_old;
    k;
    for (int i=0; i<n_dim; ++i) {
      r2[i] = r1[i] + u1[i] * segment_length_;
    }
    dl = old_segment_length - segment_length_;
    for (int i_site=2; i_site<n_sites_-1; ++i_site) {
      r2 = sites[i_site].GetPosition();
      r2_old = sites[i_site].GetPrevPosition();
      r1 = sites[i_site-1].GetPrevPosition();
      u1 = sites[i_site-1].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-2]) + SQR(segment_length_) - SQR(dl);
      k = (k>0 ? k : 0);
      k = -cos_thetas_[i_site-2]*dl + sqrt(k);
      dl = 0;
      for (int i=0; i<n_dim; ++i) {
        r2[i] = r1[i] + k * u1[i];
        dl += SQR(r2_old[i]-r2[i]);
      }
      dl = sqrt(dl);
    }
    r2 = sites[n_sites_-1].GetPosition();
    r2_old = sites[n_sites_-2].GetPrevPosition();
    r1 = sites[n_sites_-2].GetPosition();
    u1 = sites[n_sites_-3].GetOrientation();
    for (int i=0; i<n_dim; ++i) 
      r2[i] = r1[i] + segment_length_ * u1[i];
    error = old_segment_length - k - segment_length_;
    dl = ABS(error);
    j++;
    if (error!=error || j > 100) {
      std::cout << j << std::endl;
      std::cout << error << std::endl;
      error_exit("Error while adding site\n");
    }
  } while (dl > 1e-3);
  sites[n_sites_-1].SetScaledPosition(sites[n_sites_-1].GetPosition());
  length_ = segment_length_ * n_segments_;
  for (site_iterator i_site=sites.begin();
      i_site != sites.end();
      ++i_site)
    i_site->SetLength(segment_length_);
  UpdateSiteOrientations();
  CalculateAngles(false);
}

void Filament::RemoveSite() {
  int n_dim = params_->n_dim;
  UpdatePrevPositions();
  n_sites_--;
  n_segments_--;
  double old_segment_length = segment_length_;
  segment_length_ = length_/n_segments_;
  double error = 0;
  double dl = 1;
  int j=0;
  do {
    segment_length_ = segment_length_ + error/n_segments_;
    double *r, *r_old, *u;
    double k;
    dl = old_segment_length;
    for (int i_site=1; i_site<n_sites_; ++i_site) {
      r = sites[i_site].GetPosition();
      r_old = sites[i_site].GetPrevPosition();
      u = sites[i_site].GetOrientation();
      k = SQR(dl*cos_thetas_[i_site-1]) + SQR(segment_length_) - SQR(dl);
      k = (k>0 ? k : 0);
      k = -cos_thetas_[i_site-1]*dl + sqrt(k);
      for (int i=0; i<n_dim; ++i)
        r[i] = r_old[i] + k * u[i];
      dl = old_segment_length - k;
    }
    error = old_segment_length - k;
    dl = ABS(error);
    j++;
    if (error!=error || j > 100) {
      std::cout << j << std::endl;
      std::cout << error << std::endl;
      error_exit("Error while removing site\n");
    }
  } while (dl > 1e-3);
  sites.pop_back();
  length_ = segment_length_ * n_segments_;
  for (site_iterator i_site=sites.begin();
      i_site != sites.end();
      ++i_site)
    i_site->SetLength(segment_length_);
  UpdateSiteOrientations();
  CalculateAngles(false);
}

void Filament::UpdatePolyState() {
  double delta = params_->delta;
  double p_stg = params_->f_shrink_to_grow * delta;
  double p_stp = params_->f_shrink_to_pause * delta;
  double p_gts = params_->f_grow_to_shrink * delta;
  double p_gtp = params_->f_grow_to_pause * delta;
  double p_ptg = params_->f_pause_to_grow * delta;
  double p_pts = params_->f_pause_to_shrink * delta;
  double roll = gsl_rng_uniform_pos(rng_.r);
  double p_norm;
  // Modify catastrophe probabilities if the end of the filament is under a load
  if (params_->force_induced_catastrophe_flag) {
    if (wall_force_ > 0) {
      double p_factor = exp(0.0828*wall_force_);
      p_gts = (p_gts+p_gtp)*p_factor;
      p_pts = p_pts*p_factor;
    }
  }
  // Filament shrinking
  if (poly_state_ == SHRINK) {
    p_norm = p_stg + p_stp;
    if (p_norm > 1.0) 
      poly_state_ = (roll < p_stg/p_norm ? GROW : PAUSE);
    else {
      if (roll < p_stg) 
        poly_state_ = GROW;
      else if (roll < (p_stg + p_stp)) 
        poly_state_ = PAUSE;
    }
  }
  // Filament growing
  else if (poly_state_ == GROW) {
    p_norm = p_gts + p_gtp;
    if (p_norm > 1.0)
      poly_state_ = (roll < p_gts/p_norm ? SHRINK : PAUSE);
    else {
      if (roll < p_gts) 
        poly_state_ = SHRINK;
      else if (roll < (p_gts + p_gtp)) 
        poly_state_ = PAUSE;
    }
  }
  // Filament paused
  else if (poly_state_ == PAUSE) {
    p_norm = p_ptg + p_pts;
    if (p_norm > 1) 
      poly_state_ = (roll < p_ptg/p_norm ? GROW : SHRINK);
    else {
      if (roll < p_ptg) 
        poly_state_ = GROW;
      else if (roll < (p_ptg + p_pts)) 
        poly_state_ = SHRINK;
    }
  }
  // Check to make sure the filament lengths stay in the correct ranges
  if (length_ < params_->min_length)
    poly_state_ = GROW;
  else if (length_ > params_->max_length)
    poly_state_ = SHRINK;
}

void Filament::ZeroForces() {
  double fzero[3] = {0, 0, 0};
  for (int i_site=0; i_site < n_sites_; ++i_site)
    sites[i_site].SetForce(fzero);
  if (attach_state_ == ATTACHED) 
    for (int i=0; i<params_->n_dim; ++i) 
      anchor_.f_anchor[i] = 0.0;
  wall_force_ = 0.0;
}

int Filament::NSegments() {
  return n_segments_;
}

int Filament::NSites() {
  return n_sites_;
}

int Filament::GetAttachState() {
  return (int) attach_state_;
}

double *Filament::GetSitePosition(int i_site) {
  return sites[i_site].GetPosition();
}

double Filament::GetSegmentLength() {
  return segment_length_;
}

double Filament::GetDiameter() {
  return diameter_;
}

anchor_struct *Filament::GetAnchor() {
  return &anchor_;
}

double Filament::GetFrictionPerp() {
  return friction_perp_;
}

void Filament::ClearInteractions() {
  for (site_iterator i_site = sites.begin(); i_site != sites.end(); ++i_site)
    i_site->SetInteraction(false);
}

double *Filament::GetCosThetas() {
  return cos_thetas_;
}

double *Filament::GetAnchorPosition() {
  return anchor_.r_anchor;
}

double *Filament::GetAnchorScaledPosition() {
  return anchor_.s_anchor;
}

double *Filament::GetAnchorOrientation() {
  return anchor_.u_anchor;
}

void Filament::SetWallForce(double wall_force) {
  wall_force_ = wall_force;
}

double *Filament::GetAnchorForce() {
  return anchor_.f_anchor;
}
