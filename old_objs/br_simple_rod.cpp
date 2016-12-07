
#include "br_simple_rod.h"

void BrSimpleRod::Init() {
  Simple::Init();
  SetDiffusion();
  std::fill(body_frame_, body_frame_+6, 0.0);
}

void BrSimpleRod::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  Integrate();
  UpdatePeriodic();
  ClearInteractions();
}

void BrSimpleRod::UpdatePositionMP() {
  Integrate();
  UpdatePeriodic();
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
void BrSimpleRod::Integrate() {
  //Explicit calculation of Xi.F_s
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      position_[i] += 
        gamma_par_*orientation_[i]*orientation_[j]*force_[j]*delta_;
    }
    position_[i] += force_[i]*gamma_perp_*delta_;
  }
  //Add the random displacement dr(t)
  AddRandomDisplacement();
  // Update the orientation due to torques and random rotation
  UpdateOrientation();
  // Update trajectory for neighbor lists
  for (int i=0; i<n_dim_; ++i)
    dr_tot_[i] += position_[i] - prev_position_[i];
}
/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void BrSimpleRod::AddRandomDisplacement() {
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_par_);
  for (int i=0; i<n_dim_; ++i)
    position_[i] += mag * orientation_[i];
  // Then the perpendicular component(s)
  for (int j=0; j<n_dim_-1; ++j) {
    mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_perp_);
    for (int i=0; i<n_dim_; ++i)
      position_[i] += mag * body_frame_[n_dim_*j+i];
  }
  // Handle the random orientation update after updating orientation from
  // interaction torques
}

/* The orientation update is also from Yu-Guo Tao,

   u(t+dt) = u(t) + gamma_rot^-1 * T_s(t) x u(t) * dt + du(t)

   where similar to above, du(t) is the reorientation due to
   random forces, and is treated as random displacement vector(s)
   orthogonal to u(t) with std dev sqrt(2*kT*dt/gamma_rot) */
void BrSimpleRod::UpdateOrientation() {
  // First handle reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du,3);
  for (int i=0; i<n_dim_; ++i)
    orientation_[i] += du[i]*delta_/gamma_rot_;
  // Now handle the random orientation update
  for (int j=0; j<n_dim_-1; ++j) {
    double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
    for (int i=0; i<n_dim_; ++i)
      orientation_[i] += mag * body_frame_[n_dim_*j+i];
  }
  normalize_vector(orientation_, n_dim_);
}

/* calculates vector(s) orthogonal to orientation of rod */
void BrSimpleRod::GetBodyFrame() {
  if (n_dim_==2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  }
  else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]),3);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]),3);
    normalize_vector(&(body_frame_[0]),n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]),3);
  }
}

/* Initialize diffusion coefficients and std dev for random numbers */
void BrSimpleRod::SetDiffusion() {
  double logLD = log(length_/diameter_);
  gamma_par_ = 2.0*length_ / (3.0*logLD);
  gamma_perp_ = 2.0*gamma_par_;
  gamma_rot_ = length_*length_*length_ / (9.0*logLD);
  rand_sigma_par_ = sqrt(2.0*delta_/gamma_par_);
  rand_sigma_perp_ = sqrt(2.0*delta_/gamma_perp_);
  rand_sigma_rot_ = sqrt(2.0*delta_/gamma_rot_);
}

//void BrSimpleRodSpecies::InitPotentials (system_parameters *params) {
  //AddPotential(SID::br_simple_rod, SID::br_simple_rod, 
      //// Set br_rod-br_rod interaction
      //new WCA(params->lj_epsilon,params->rod_diameter,
        //space_, pow(2, 1.0/6.0)*params->rod_diameter));
//}
