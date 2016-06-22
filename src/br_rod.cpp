#include "br_rod.h"

void BrRod::Init() {
  InsertRandom(0.5*length_+diameter_);
  printf("Inserted rod: \n");
  printf("  pos : {%f, %f, %f}\n",position_[0],position_[1],position_[2]);
  printf("  u : {%f, %f, %f}\n",orientation_[0],orientation_[1],orientation_[2]);
  printf("  l : %f\n",length_);
  printf("  d : %f\n",diameter_);
  SetDiffusion();
  std::fill(body_frame_, body_frame_+6, 0.0);
  // Init bond lengths and diameter
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetLength(child_length_);
    bond->SetDiameter(diameter_);
  }
  // Set positions for sites and bonds
  UpdateSiteBondPositions();
}

void BrRod::UpdatePositionMP() {
  Integrate();
  UpdatePeriodic();
  // Update end site positions for tracking trajectory for neighbors
  UpdateSiteBondPositions();
}

void BrRod::UpdateSiteBondPositions() {
  // First set prev positions for sites
  elements_[0].SetPrevPosition(elements_[0].GetPosition());
  elements_[1].SetPrevPosition(elements_[1].GetPosition());
  double pos[3];
  // then update site positions based on new COM and orientation
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] - 0.5 * length_ * orientation_[i];
  elements_[0].SetPosition(pos);
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] + 0.5 * length_ * orientation_[i];
  elements_[1].SetPosition(pos);
  // update trajectories of end sites for neighbor lists
  elements_[0].AddDr();
  elements_[1].AddDr();
  // Set site of first bond COM and update remaining COMs
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] + 0.5*(child_length_-length_)*orientation_[i];
  // XXX Get rid of prints, k, JMM
  int k=1;
  double u[3];
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetPosition(pos);
    if (k%2 == 0)
      bond->SetOrientation(orientation_);
    else {
      for (int i=0; i<n_dim_; ++i)
        u[i] = -orientation_[i];
      bond->SetOrientation(u);
    }
    // Set next bond COM
    printf("bond %d:\n",k);
    printf("  pos : {%f, %f, %f}\n",pos[0],pos[1],pos[2]);
    printf("  u : {%f, %f, %f}\n",orientation_[0],orientation_[1],orientation_[2]);
    printf("  l : %f\n",bond->GetLength());
    printf("  d : %f\n",bond->GetDiameter());
    for (int i=0; i<n_dim_; ++i)
      pos[i] += orientation_[i] * child_length_;
    k++;
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
void BrRod::Integrate() {
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
}

double const * const BrRod::GetDrTot() {
  double dr_max=0;
  for (auto site=elements_.begin(); site!= elements_.end(); ++site) {
    double dr_mag = 0;
    double const * const dr = site->GetDrTot();
    for (int i=0; i<n_dim_; ++i)
      dr_mag += dr[i]*dr[i];
    if (dr_mag > dr_max) {
      dr_max = dr_mag;
      std::copy(dr, dr+3, dr_tot_);
    }
  }
  return dr_tot_;
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void BrRod::AddRandomDisplacement() {
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
void BrRod::UpdateOrientation() {
  // First handle reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, n_dim_);
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
void BrRod::GetBodyFrame() {
  if (n_dim_==2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  }
  else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]),n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

/* Initialize diffusion coefficients and std dev for random numbers */
void BrRod::SetDiffusion() {
  double logLD = log(length_/diameter_);
  gamma_par_ = 2.0*length_ / (3.0*logLD);
  gamma_perp_ = 2.0*gamma_par_;
  gamma_rot_ = length_*length_*length_ / (9.0*logLD);
  rand_sigma_par_ = sqrt(2.0*delta_/gamma_par_);
  rand_sigma_perp_ = sqrt(2.0*delta_/gamma_perp_);
  rand_sigma_rot_ = sqrt(2.0*delta_/gamma_rot_);
}

void BrRodSpecies::InitPotentials (system_parameters *params) {
  AddPotential(SID::br_simple_rod, SID::br_simple_rod, 
      // Set br_rod-br_rod interaction
      new WCA(params->lj_epsilon,params->rod_diameter,
        space_, pow(2, 1.0/6.0)*params->rod_diameter));
}

void BrRod::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)
    bond->Draw(graph_array);
}

