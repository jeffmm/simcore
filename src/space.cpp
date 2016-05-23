#include "space.h"

// Constructor for array reference construction
SpaceProperties::SpaceProperties() {}

void SpaceProperties::Init(system_parameters *params, long seed) {

  params_ = params;
  n_dim_ = params_->n_dim;
  n_periodic_ = params_->n_periodic;
  radius_ = params_->system_radius;
  d_radius_ = 0;
  m_d_dist_ = 0;
  r_cutoff_ = params_->r_cutoff_boundary;
  // Make sure space_type is recognized
  if (radius_ <= 0)
    error_exit("ERROR: Mother cell radius is not a positive number!\n");
  if (params_->boundary_type >= 0 && params_->boundary_type < 2) {
    boundary_type_ = (params_->boundary_type == 0 ? SPHERE : BOX);
  }
  else if (params_->boundary_type == 2) {
    d_radius_ = params_->daughter_radius;
    m_d_dist_ = params_->mother_daughter_dist;
    if (d_radius_ > radius_) 
      error_exit("ERROR: Daughter cell radius larger than mother cell radius!\n");
    if (d_radius_ <= 0 || m_d_dist_ <= 0)
      error_exit("ERROR: Snowman boundary type selected but either daughter cell radius\n\
          or mother/daughter cell separation is not a positive number!\n");
    boundary_type_ = SNOWMAN;
  } 
  else {
    error_exit("ERROR: Invalid boundary type given.\n");
  }
  InitUnitCell();
  CalculateVolume();
  rng_.init(seed);
  UpdateSpaceStruct();
}

void SpaceProperties::Clear() {
  rng_.clear();
  ClearUnitCell();
}


void SpaceProperties::InitUnitCell() {
  // This needs to be changed to account for periodic boundary conditions
  uc_ = new double*[n_dim_];
  uc_inv_ = new double*[n_dim_];
  a_ = new double*[n_dim_];
  b_ = new double*[n_dim_];
  a_perp_ = new double[n_dim_];

  for(int i=0; i<n_dim_; ++i) {
      uc_[i] = new double[n_dim_]; 
      uc_inv_[i] = new double[n_dim_];
      a_[i] = new double[n_dim_];
      b_[i] = new double[n_dim_];
  }

  double h_param1 = 2.0 * radius_;
  double h_param2 = m_d_dist_ + d_radius_ + radius_;
  double h_major = (h_param1 > h_param2 ? h_param1 : h_param2);
  double h_minor;
  if (boundary_type_ == SNOWMAN)
    h_minor = (h_param1 > h_param2 ? h_param2 : h_param1);
  else
    h_minor = h_major;
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      uc_[i][j] = (i==j ? 1 : 0) * h_minor;
    }
  }
  uc_[n_dim_-1][n_dim_-1] = h_major;

  /* Compute inverse unit cell matrix. */
  if (n_dim_==2)
    invert_sym_2d_matrix(uc_, uc_inv_);
  if (n_dim_==3)
    invert_sym_3d_matrix(uc_, uc_inv_);

  /* Compute unit cell volume. */
  uc_volume_ = determinant(n_dim_, uc_);

  /* Set up direct and reciprocal lattice vectors. */
  for (int i = 0; i < n_dim_; ++i)
    for (int j = 0; j < n_dim_; ++j) {
      a_[i][j] = uc_[j][i];
      b_[i][j] = uc_inv_[i][j];
    }

  /* Compute perpendicular distances between opposite unit cell faces. */
  for (int i = 0; i < n_dim_; ++i) {
    double b_mag2 = 0.0;
    for (int j = 0; j < n_dim_; ++j)
      b_mag2 += SQR(b_[i][j]);
    double b_mag = sqrt(b_mag2);
    a_perp_[i] = 1.0 / b_mag;
  }
}

void SpaceProperties::ClearUnitCell() {
  for (int i=0; i<n_dim_; ++i) {
    delete[] uc_[i];
    delete[] uc_inv_[i];
    delete[] a_[i];
    delete[] b_[i];
  }
  delete[] uc_;
  delete[] uc_inv_;
  delete[] a_perp_;
  delete[] a_;
  delete[] b_;
}

void SpaceProperties::CalculateVolume() {

  if (boundary_type_ == SPHERE) {
    v_ratio_ = 0;
    intersect_height_ = 0;
    intersect_radius_ = 0;
    if (n_dim_ == 2)
      volume_ = M_PI * SQR(radius_);
    else
      volume_ = 4.0/3.0*M_PI*CUBE(radius_);
  }
  else if (boundary_type_ == BOX) {
    v_ratio_ = 0;
    intersect_height_ = 0;
    intersect_radius_ = 0;
    //FIXME This is not correct for cell lists with periodic boundary. Should not hurt calculations though, just not efficient.
    if (n_dim_ == 2)
      volume_ = SQR(2*radius_);
    else
      volume_ = CUBE(2*radius_);
  }
  else if (boundary_type_ == SNOWMAN) {
    double R = radius_;
    double r = d_radius_;
    double d = m_d_dist_;
    if (n_dim_ == 2) {
      volume_ = M_PI * SQR(R) + M_PI * SQR(r)
               - SQR(r) * acos((SQR(d)+SQR(r)-SQR(R))/(2*d*r))
               - SQR(R) * acos((SQR(d)+SQR(R)-SQR(r))/(2*d*R))
               + 0.5 * sqrt((R+r-d)*(r+d-R)*(R+d-r)*(R+r+d));
      v_ratio_ = SQR(r) / (SQR(r) + SQR(R));
    }
    else {
      volume_ = 4.0/3.0*M_PI*CUBE(R) + 4.0/3.0*M_PI*CUBE(r)
                - M_PI * SQR(R+r-d)
                  * (SQR(d) + 2*d*r - 3*SQR(r) + 2*d*R + 6*r*R - 3*SQR(R))
                  / (12*d);
      v_ratio_ = CUBE(r) / (CUBE(R) + CUBE(r));
    }
    intersect_height_ = ( SQR(R)+SQR(d)-SQR(r) )/(2*d);
    intersect_radius_ = R * sqrt(1 - SQR( (SQR(R)+SQR(d)-SQR(r))/(2*d*R) ));
  }
}

void SpaceProperties::RandomCoordinate(double *vec) {
  RandomCoordinate(vec, 0);
}

void SpaceProperties::RandomCoordinate(double *vec, double buffer) {
  double mag;
  if (radius_ - r_cutoff_ - buffer < 0) 
    error_exit("ERROR: Cannot return random coordinate in system volume with given buffer\n");
  if (boundary_type_ == SPHERE) {
    generate_random_unit_vector(n_dim_, vec, rng_.r);
    mag = gsl_rng_uniform_pos(rng_.r) * (radius_ - r_cutoff_ - buffer);
    for (int i=0; i<n_dim_; ++i) {
      vec[i] *= mag;
    }
  }
  else if (boundary_type_ == BOX) {
    for (int i=0; i<n_dim_; ++i)
      vec[i] = (2.0*gsl_rng_uniform_pos(rng_.r)-1.0) * (radius_ - r_cutoff_ - buffer);
  }
  else if (boundary_type_ == SNOWMAN) {
    double roll = gsl_rng_uniform_pos(rng_.r);
    mag = gsl_rng_uniform_pos(rng_.r);
    generate_random_unit_vector(n_dim_, vec, rng_.r);
    if (roll < v_ratio_) {
      // Place coordinate in daughter cell
      mag *= (d_radius_ - r_cutoff_ - buffer);
      for (int i=0; i<n_dim_; ++i) {
        vec[i] *= mag;
      }
      vec[n_dim_-1] += m_d_dist_;
    }
    else {
      mag *= (radius_ - r_cutoff_ - buffer);
      for (int i=0; i<n_dim_; ++i) {
        vec[i] *= mag;
      }
    }
  }
}

// returns true if the coordinate is within 
// the volume of the system boundary
bool SpaceProperties::CheckInBounds(double *vec, double buffer) {

  if (boundary_type_ == SPHERE) {
    double mag2 = 0.0;
    for (int i=0; i<n_dim_; ++i) 
      mag2 += SQR(vec[i]);
    if (mag2 >= SQR(radius_-r_cutoff_-buffer)) 
      return false;
  }
  else if (boundary_type_ == BOX) {
    for (int i=0; i<n_dim_; ++i) 
      if (vec[i] > (radius_-r_cutoff_-buffer) || vec[i] < -(radius_-r_cutoff_-buffer))
        return false;
  }
  else if (boundary_type_ == SNOWMAN) {
    double mag2 = 0.0;
    double dmag2 = 0.0;
    for (int i=0; i<n_dim_-1; ++i) {
      mag2 += SQR(vec[i]);
      dmag2 += SQR(vec[i]);
    }
    mag2 += SQR(vec[n_dim_-1]);
    dmag2 += SQR(vec[n_dim_-1] - m_d_dist_);
    if (mag2 > SQR(radius_-r_cutoff_-buffer) && dmag2 > SQR(d_radius_-r_cutoff_-buffer)) {
      return false;
    }
  }
  return true;
}

// Returns true if the segment defined by the line between 
// the points vec1 and vec2 falls outside the snowman boundary
bool SpaceProperties::CheckSegmentInBounds(double *vec1, double *vec2, double buffer) {

  // If both endpoints are on the same side of the mother/daughter cell
  // intersection, the segment must be in bounds if the sites are in bounds
  double *u_vec = new double[n_dim_];
  if ( SIGNOF(vec1[n_dim_-1] - intersect_height_) != SIGNOF( vec2[n_dim_-1] - intersect_height_) ) {
    // Check for point that passes through the plane at mother, daughter cell intersection
    double vec_mag = 0.0;
    for (int i=0; i<n_dim_; ++i) {
      vec_mag += SQR(vec1[i] - vec2[i]);
    }
    for (int i=0; i<n_dim_; ++i) 
      u_vec[i] = (vec1[i] - vec2[i])/vec_mag; 
    // For the (improbable) case that the filament is nearly parallel to the plane
    if ( ABS(u_vec[n_dim_-1]) < 1.0e-8 ) {
      delete[] u_vec;
      return true;
    }
    double d = (intersect_height_ - vec1[n_dim_-1])/u_vec[n_dim_-1];
    double intersection_mag = 0.0;
    for (int i=0; i<n_dim_-1; ++i) {
      intersection_mag += SQR(vec1[i] + d * u_vec[i]);
    }
    if (intersection_mag > SQR(intersect_radius_ - r_cutoff_ - buffer)) {
      delete[] u_vec;
      return false;
    }
  }
  delete[] u_vec;
  return true;
}

int const SpaceProperties::GetDim() {
  return n_dim_;
}

int const SpaceProperties::GetPeriodic() {
  return n_periodic_;
}

double const SpaceProperties::GetRadius() {
  return radius_;
}

double const SpaceProperties::GetDRadius() {
  return d_radius_;
}

double const SpaceProperties::GetMDDist() {
  return m_d_dist_;
}

double const SpaceProperties::GetVolume() {
  return volume_;
}

double const * const * const SpaceProperties::GetUnitCell() {
  return uc_;
}

double const * const * const SpaceProperties::GetUnitCellInv(){
    return uc_inv_;
}

boundary_type_t SpaceProperties::GetType() {
  return boundary_type_;
}

std::string SpaceProperties::GetTypeString() {
  std::string type;
  switch (boundary_type_) {
    case 0:
      type = "sphere";
      break;
    case 1:
      type = "cube";
      break;
    case 2:
      type = "snowman";
      break;
    default:
      type = "none";
      break;
  }
  return type;
}

double const * const SpaceProperties::GetAPerp() {
    return a_perp_;
}
double const * const * const SpaceProperties::GetA() {
    return a_;
}
double const * const * const SpaceProperties::GetB() {
    return b_;
}

double const SpaceProperties::GetIntersectHeight() {
  return intersect_height_;
}

double const SpaceProperties::GetIntersectRadius() {
  return intersect_radius_;
}

void SpaceProperties::UpdateSpaceStruct() {
  s_struct.n_dim = n_dim_;
  s_struct.n_periodic = n_periodic_;
  s_struct.type = GetTypeString();
  if (GetType() == 2) {
    s_struct.bud = true;
    s_struct.bud_height = m_d_dist_;
    s_struct.bud_radius = d_radius_;
  }
  else {
    s_struct.bud = false;
  }
  s_struct.radius = radius_;
  s_struct.unit_cell = uc_;
  s_struct.unit_cell_inv = uc_inv_;
  s_struct.a = a_;
  s_struct.b = b_;
  s_struct.a_perp = a_perp_;
}

space_struct * SpaceProperties::GetStruct() {
  return &s_struct;
}

//void SpaceProperties::InsertRandom(double pos[3], double buffer) {
  //double *p = new double[3];
  //RandomCoordinate(p,buffer);
  //memcpy(pos,p,3*sizeof(double));
  //delete p;
//}

