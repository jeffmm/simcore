#include "space.h"

// Constructor for array reference construction
Space::Space() {}

void Space::Init(system_parameters *params) {

  params_ = params;
  n_dim_ = params_->n_dim;
  n_periodic_ = params_->n_periodic;
  delta_ = params_->delta;
  compressibility_ = params_->compressibility;
  pressure_time_ = params_->pressure_time;
  // Make sure that n_periodic <= n_dim
  if (n_periodic_ > n_dim_) 
    n_periodic_ = params_->n_periodic = n_dim_;
  radius_ = params_->system_radius;
  if (radius_ <= 0)
    error_exit("ERROR: System radius is not a positive number!\n");
  bud_radius_ = 0;
  bud_height_ = 0;
  //r_cutoff_ = params_->r_cutoff_boundary;
  constant_pressure_ = (params_->constant_pressure == 1 ? true : false);
  constant_volume_ = (params_->constant_volume == 1 ? true : false);
  target_pressure_ = params_->target_pressure;
  target_radius_ = params_->target_radius;
  update_ = false;

  // Make sure space_type is recognized
  switch (params_->boundary_type) {
    case 0:
      boundary_type_ = BOX;
      boundary_type_string_ = "box";
      break;
    case 1:
      boundary_type_ = SPHERE;
      boundary_type_string_ = "sphere";
      break;
    case 2:
      // No periodicity allowed in budding yeast boundary type
      n_periodic_ = params_->n_periodic = 0;
      // Get bud parameters
      bud_radius_ = params_->bud_radius;
      bud_height_ = params_->bud_height;
      if (bud_radius_ > radius_) 
        error_exit("ERROR: Bud radius larger than parent cell radius!\n");
      if (bud_radius_ <= 0 || bud_height_ <= 0)
        error_exit("ERROR: Budding yeast boundary type selected but either bud radius\n\
            or bud height is not a positive number!\n");
      boundary_type_ = BUDDING;
      boundary_type_string_ = "budding";
    default:
      error_exit("ERROR: Boundary type %d not recognized!\n",params_->boundary_type);
  }
  InitUnitCell();
  CalculateVolume();
  InitSpaceStruct();
}

void Space::InitUnitCell() {

  double h_param1 = 2.0 * radius_;
  double h_major = h_param1;
  double h_minor = h_major;
  if (boundary_type_ == BUDDING) {
    double h_param2 = bud_height_ + bud_radius_ + radius_;
    h_major = (h_param1 > h_param2 ? h_param1 : h_param2);
    h_minor = (h_param1 > h_param2 ? h_param2 : h_param1);
  }
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      unit_cell_[i*n_dim_+j] = (i==j ? 1 : 0) * h_minor;
    }
  }
  unit_cell_[n_dim_*n_dim_-1] = h_major;
  // Compute unit cell quantities
  CalculateUnitCellQuantities();

}

void Space::CalculateUnitCellQuantities() {
  // Calculate unit cell volume, which is determinant of unit cell
  double determinant;
  if (n_dim_==2) {
    determinant = unit_cell_[0]*unit_cell_[3] - unit_cell_[2]*unit_cell_[1];
  }
  else if (n_dim_==3) {
    determinant = unit_cell_[0]*((unit_cell_[4]*unit_cell_[8]) - (unit_cell_[7]*unit_cell_[5])) -unit_cell_[1]*(unit_cell_[3]*unit_cell_[8] - unit_cell_[6]*unit_cell_[5]) + unit_cell_[2]*(unit_cell_[3]*unit_cell_[7] - unit_cell_[6]*unit_cell_[4]);
  }
  unit_cell_volume_ = determinant;
  /* Compute inverse unit cell matrix. */
  if (n_dim_==2)
    invert_sym_2d_matrix(unit_cell_, unit_cell_inv_);
  if (n_dim_==3)
    invert_sym_3d_matrix(unit_cell_, unit_cell_inv_);
  /* Compute unit cell volume, which is determinant of uc matrix */

  /* Set up direct and reciprocal lattice vectors. */
  for (int i = 0; i < n_dim_; ++i)
    for (int j = 0; j < n_dim_; ++j) {
      a_[n_dim_*i+j] = unit_cell_[n_dim_*j+i];
      b_[n_dim_*i+j] = unit_cell_inv_[n_dim_*i+j];
    }

  /* Compute perpendicular distances between opposite unit cell faces. */
  for (int i = 0; i < n_dim_; ++i) {
    double b_mag2 = 0.0;
    for (int j = 0; j < n_dim_; ++j)
      b_mag2 += SQR(b_[n_dim_*i+j]);
    double b_mag = sqrt(b_mag2);
    a_perp_[i] = 1.0 / b_mag;
  }
}


void Space::UpdateSpace() {
  // No space update for budding yeast yet
  if (boundary_type_ == BUDDING)
    return;
  UpdateUnitCell();
  UpdateVolume();
  UpdateSpaceStruct();
}

/* XXX Only works for periodic boundary conditions so far, need to add volume
       updating for fixed boundaries (like enclosed boxes or spheres) */
void Space::ConstantPressure() {
  update_ = true;
  pressure_ = s_struct.pressure;
  // If target pressure is zero, let first pressure calculation set target pressure
  if (target_pressure_ == 0)
    target_pressure_ = pressure_; 
  // If target pressure is vastly different from current pressure, update scaling matrix
  if (ABS(pressure_ - target_pressure_) > 1e-4) 
    CalculateScalingMatrix();
  else
    update_ = false;
  printf("pressure: %2.8f\ntarget_pressure: %2.8f\nradius: %2.8f\n",pressure_,target_pressure_,radius_);
}

void Space::ConstantVolume() {
  // If target radius is <= zero, disable change of system
  if (target_radius_ <= 0)
    target_radius_ = radius_;
  // If target pressure is vastly different from current pressure, update scaling matrix
  if (ABS(target_radius_ - radius_) > 2*delta_) {
    update_ = true;
    CalculateScalingMatrix();
  }
  else
    update_ = false;
  printf("target_radius: %2.8f\nradius: %2.8f\n",target_radius_,radius_);
}

// Use scaling matrix to scale unit cell
void Space::UpdateUnitCell() {
  std::copy(unit_cell_, unit_cell_+9, prev_unit_cell_);
  std::fill(unit_cell_, unit_cell_+9, 0);
  for (int i=0; i<n_dim_; ++i)
    for (int j=0; j<n_dim_; ++j)
      for (int k=0; k<n_dim_; ++k)
        unit_cell_[n_dim_*i+j] += mu_[n_dim_*i+k] * prev_unit_cell_[n_dim_*k+j];
  CalculateUnitCellQuantities();
}

/* Warning: Only does orthogonal scaling 
   (doesn't work for unit cells with non-zero off-diagonal elements) */
void Space::CalculateScalingMatrix() {
  double time_const = compressibility_*delta_/(n_dim_*pressure_time_);
  std::copy(s_struct.pressure_tensor,s_struct.pressure_tensor+9,pressure_tensor_);
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      if (i==j) {
        if (constant_pressure_)
          mu_[n_dim_*i+j] = 1.0 - time_const*(target_pressure_ - pressure_tensor_[n_dim_*i+j]);
        else if (constant_volume_)
          mu_[n_dim_*i+j] = 1.0 - time_const*(0.5*unit_cell_[n_dim_*i+j] - target_radius_);
      }
      else 
        mu_[n_dim_*i+j] = 0.0; // no non-ortho scaling
    }
  }
}

void Space::CalculateRadius() {
  if (boundary_type_ == BOX) {
    if (n_dim_ == 3)
      radius_ = 0.5*pow(volume_, 1.0/3.0);
    else if (n_dim_ == 2)
      radius_ = 0.5*sqrt(volume_);
  }
  else if (boundary_type_ == SPHERE) {
    if (n_dim_ == 3)
      radius_ = pow(3.0*volume_/(4.0*M_PI), 1.0/3.0);
    else if (n_dim_ == 2)
      radius_ = sqrt(volume_/M_PI);
  }
}

void Space::UpdateVolume() {
  if (boundary_type_ == BOX) {
    volume_ = unit_cell_volume_;
    CalculateRadius();
  }
  else if (boundary_type_ == SPHERE) {
    error_exit("ERROR: Updating volume for spherical boundary types not implemented.\n");
  }
  else if (boundary_type_ == BUDDING) {
    error_exit("ERROR: Updating volume for budding yeast boundary types not implemented.\n");
  }
}

void Space::CalculateVolume() {
  if (boundary_type_ == BOX) {
    v_ratio_ = 0;
    neck_height_ = 0;
    neck_radius_ = 0;
    //FIXME This is not correct for cell lists with periodic boundary. Should not hurt calculations though, just not efficient.
    if (n_dim_ == 2)
      volume_ = SQR(2*radius_);
    else
      volume_ = CUBE(2*radius_);
  }
  else if (boundary_type_ == SPHERE) {
    v_ratio_ = 0;
    neck_height_ = 0;
    neck_radius_ = 0;
    if (n_dim_ == 2)
      volume_ = M_PI * SQR(radius_);
    else
      volume_ = 4.0/3.0*M_PI*CUBE(radius_);
  }
  else if (boundary_type_ == BUDDING) {
    double R = radius_;
    double r = bud_radius_;
    double d = bud_height_;
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
    neck_height_ = ( SQR(R)+SQR(d)-SQR(r) )/(2*d);
    neck_radius_ = R * sqrt(1 - SQR( (SQR(R)+SQR(d)-SQR(r))/(2*d*R) ));
  }
}

void Space::InitSpaceStruct() {
  s_struct.n_dim = n_dim_;
  s_struct.n_periodic = n_periodic_;
  s_struct.type = boundary_type_string_;
  if (boundary_type_ == BUDDING) {
    s_struct.bud = true;
    s_struct.bud_height = bud_height_;
    s_struct.bud_radius = bud_radius_;
  }
  else {
    s_struct.bud = false;
  }
  s_struct.unit_cell = unit_cell_;
  s_struct.unit_cell_inv = unit_cell_inv_;
  s_struct.a = a_;
  s_struct.b = b_;
  s_struct.a_perp = a_perp_;
  s_struct.mu = mu_;
  UpdateSpaceStruct();
}

// Only need to update changing variables (but not pointers)
void Space::UpdateSpaceStruct() {
  s_struct.radius = radius_;
  s_struct.volume = volume_;
}

space_struct * Space::GetStruct() {
  return &s_struct;
}

