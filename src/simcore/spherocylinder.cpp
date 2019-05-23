#include "spherocylinder.h"

Spherocylinder::Spherocylinder() : Object() {
  color_ = params_->spherocylinder.color;
  draw_ = draw_type::_from_string(params_->spherocylinder.draw_type.c_str());
  diameter_ = params_->spherocylinder.diameter;
  length_ = params_->spherocylinder.length;
  is_midstep_ = params_->spherocylinder.midstep;
  std::fill(body_frame_, body_frame_+6, 0.0);
  SetDiffusion();
}

void Spherocylinder::Init() {
  do {
    InsertSpherocylinder();
  } while (CheckBounds());
}

void Spherocylinder::InsertSpherocylinder() {
  if (params_->spherocylinder.insertion_type.compare("random") == 0) {
    InsertRandom();
  }
  else if (params_->spherocylinder.insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
  }
  else if (params_->spherocylinder.insertion_type.compare("centered_random") == 0) {
    std::fill(position_,position_+3,0.0);
    generate_random_unit_vector(n_dim_,orientation_,rng_.r);
  }
  else if (params_->spherocylinder.insertion_type.compare("centered_oriented") == 0 ) {
    std::fill(position_,position_+3,0.0);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
  }
  else {
    error_exit("Spherocylinder insertion type not recognized!");
  }
}

void Spherocylinder::UpdatePosition() {
  SetPrevPosition(position_);
  ApplyForcesTorques();
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
void Spherocylinder::Integrate() {
  double delta = (is_midstep_ ? 0.5*delta_ : delta_);
  //Explicit calculation of Xi.F_s
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      position_[i] +=
        gamma_par_*orientation_[i]*orientation_[j]*force_[j]*delta;
    }
    position_[i] += force_[i]*gamma_perp_*delta;
  }
  // Reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, 3); // ndim=3 since torques
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] += du[i]*delta/gamma_rot_;
  }
  //Add the random displacement dr(t)
  AddRandomDisplacement();
  //Update the orientation due to torques and random rotation
  AddRandomReorientation();
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void Spherocylinder::AddRandomDisplacement() {
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  double mag = gsl_ran_gaussian_ziggurat(rng_.r, diffusion_par_);
  for (int i=0; i<n_dim_; ++i)
    position_[i] += mag * orientation_[i];
  // Then the perpendicular component(s)
  for (int j=0; j<n_dim_-1; ++j) {
    mag = gsl_ran_gaussian_ziggurat(rng_.r, diffusion_perp_);
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
void Spherocylinder::AddRandomReorientation() {
  // Now handle the random orientation update
  for (int j=0; j<n_dim_-1; ++j) {
    double mag = gsl_ran_gaussian_ziggurat(rng_.r, diffusion_rot_);
    for (int i=0; i<n_dim_; ++i) {
      orientation_[i] += mag * body_frame_[n_dim_*j+i];
    }
  }
  normalize_vector(orientation_, n_dim_);
}

void Spherocylinder::ApplyForcesTorques() {}

void Spherocylinder::SetDiffusion() {
  // Sets the diffusion in accordance to Lowen, Phys. Rev. E, 1994.
  double L = length_ + diameter_;
  double p = L/diameter_;
  double log_p = log(p);
  gamma_par_ = 2.0*L/3.0/(log_p - 0.207 + 0.980/p - 0.133/SQR(p));
  gamma_perp_ = 4.0*L/3.0/(log_p + 0.839 + 0.185/p + 0.233/SQR(p));
  gamma_rot_ = CUBE(L)/9.0/(log_p - 0.662 + 0.917/p - 0.050/SQR(p));
  double delta = (is_midstep_ ? 0.5*delta_ : delta_);
  diffusion_par_ = sqrt(2*delta/gamma_par_);
  diffusion_perp_ = sqrt(2*delta/gamma_perp_);
  diffusion_rot_ = sqrt(2*delta/gamma_rot_);
}

void Spherocylinder::GetBodyFrame() {
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

void SpherocylinderSpecies::InitAnalysis() {
  if (params_->spherocylinder.diffusion_analysis) {
    InitDiffusionAnalysis();
  }
}

void SpherocylinderSpecies::RunAnalysis() {
  if (params_->spherocylinder.diffusion_analysis) {
    DiffusionAnalysis();
  }
}

void SpherocylinderSpecies::FinalizeAnalysis() {
  if (params_->spherocylinder.diffusion_analysis) {
    FinalizeDiffusionAnalysis();
  }
}

void SpherocylinderSpecies::InitDiffusionAnalysis() {
  if (n_members_ == 1) {
    warning("Diffusion analysis incompatible with simulations of 1 species member. Aborting diffusion analysis");
    params_->spherocylinder.diffusion_analysis = 0;
    return;
  }
  n_samples_ = params_->spherocylinder.n_diffusion_samples;
  time_ = 0;
  int n_data = params_->n_steps / params_->spherocylinder.n_posit;
  time_avg_interval_ = n_data / n_samples_;
  if (time_avg_interval_ < 1) {
    error_exit("Something went wrong in InitDiffusionAnalysis!");
  }
  pos0_ = new double*[n_members_];
  u0_ = new double*[n_members_];
  for (int i=0;i<n_members_;++i) {
    pos0_[i] = new double[params_->n_dim];
    u0_[i] = new double[params_->n_dim];
  }
  vcf_ = new double[time_avg_interval_];
  msd_ = new double[time_avg_interval_];
  vcf_err_ = new double[time_avg_interval_];
  msd_err_ = new double[time_avg_interval_];
  std::fill(msd_,msd_+time_avg_interval_,0.0);
  std::fill(vcf_,vcf_+time_avg_interval_,0.0);
  std::fill(msd_err_,msd_err_+time_avg_interval_,0.0);
  std::fill(vcf_err_,vcf_err_+time_avg_interval_,0.0);
  UpdateInitPositions();
}


void SpherocylinderSpecies::DiffusionAnalysis() {
  // Calculate MSD
  CalculateMSD();
  CalculateVCF();
  if (++time_ == time_avg_interval_) {
    time_ = 0;
    UpdateInitPositions();
  }
}

void SpherocylinderSpecies::UpdateInitPositions() {
  for (int i=0;i<n_members_; ++i) {
    double const * const position0 = members_[i].GetPosition();
    double const * const orientation0 = members_[i].GetOrientation();
    for (int j=0;j<params_->n_dim;++j) {
      pos0_[i][j] = position0[j];
      u0_[i][j] = orientation0[j];
    }
  }
}

void SpherocylinderSpecies::CalculateMSD() {
  double avg_sqr_dist = 0.0;
  double avg_sqr_dist_sqr = 0.0;
  for (int i=0; i<n_members_; ++i) {
    double sqr_diff = 0;
    double const * const position = members_[i].GetPosition();
    for (int j=0; j<params_->n_dim; ++j) {
      double r_diff = position[j] - pos0_[i][j];
      sqr_diff += SQR(r_diff);
    }
    avg_sqr_dist += sqr_diff;
    avg_sqr_dist_sqr += SQR(sqr_diff);
  }
  avg_sqr_dist/=n_members_;
  avg_sqr_dist_sqr/=n_members_;
  double stdev2 = avg_sqr_dist_sqr - SQR(avg_sqr_dist);
  if (stdev2 < 0) {
    error_exit("Something was negative in diffusion analysis when it shouldn't have been!\n");
  }
  msd_[time_] += avg_sqr_dist/stdev2;
  msd_err_[time_] += 1.0/stdev2;
}

void SpherocylinderSpecies::CalculateVCF() {
  double avg_udotu0 = 0.0;
  double avg_udotu0_sqr = 0.0;
  for (int i=0; i<n_members_; ++i) {
    double udotu0 = 0.0;
    double const * const orientation = members_[i].GetOrientation();
    for (int j=0; j<params_->n_dim; ++j) {
      udotu0 += orientation[j] * u0_[i][j];
    }
    avg_udotu0 += udotu0;
    avg_udotu0_sqr += udotu0*udotu0;
  }
  avg_udotu0/=n_members_;
  avg_udotu0_sqr/=n_members_;
  double stdev2 = avg_udotu0_sqr - SQR(avg_udotu0);
  vcf_[time_] += avg_udotu0/stdev2;
  vcf_err_[time_] += 1.0/stdev2;
}

void SpherocylinderSpecies::FinalizeDiffusionAnalysis() {
  for (int t=0; t<time_avg_interval_; ++t) {
    msd_[t] = msd_[t]/msd_err_[t];
    vcf_[t] = vcf_[t]/vcf_err_[t];
    msd_err_[t] = 1.0/sqrt(n_members_*msd_err_[t]);
    vcf_err_[t] = 1.0/sqrt(n_members_*vcf_err_[t]);
  }
  std::string fname = params_->run_name;
  fname.append("_spherocylinder.diffusion");
  diff_file_.open(fname, std::ios::out);
  diff_file_ << "length diameter n_dim delta n_steps n_posit n_objs n_samples\n";
  diff_file_ << params_->spherocylinder.length << " " << params_->spherocylinder.diameter 
    << " " << params_->n_dim << " " << params_->delta << " " << params_->n_steps << " " 
    << params_->spherocylinder.n_posit << " " << n_members_ << " " << n_samples_ << "\n";
  diff_file_ << "time msd msd_err vcf vcf_err\n";
  diff_file_ << "0.0 0.0 0.0 1.0 0.0\n";
  diff_file_.precision(16);
  diff_file_.setf(std::ios::fixed);
  diff_file_.setf(std::ios::showpoint);
  double midterm = (midstep_ ? 0.5 : 1);
  for (int t=0; t<time_avg_interval_; ++t) {
    diff_file_ << midterm*(t+1)*params_->delta*params_->spherocylinder.n_posit << " " << msd_[t] << " " << msd_err_[t]
      << " " << vcf_[t] << " " << vcf_err_[t] << "\n";
  }
  diff_file_.close();

  for (int i=0;i<n_members_;++i) {
    delete pos0_[i];
    delete u0_[i];
  }
  delete pos0_;
  delete u0_;
}

