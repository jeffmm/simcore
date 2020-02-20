#include "simcore/spherocylinder_species.hpp"

SpherocylinderSpecies::SpherocylinderSpecies(unsigned long seed)
    : Species(seed) {
  SetSID(species_id::spherocylinder);
}
void SpherocylinderSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  midstep_ = sparams_.midstep;
}

void SpherocylinderSpecies::InitAnalysis() {
  if (sparams_.diffusion_analysis) {
    InitDiffusionAnalysis();
  }
}

void SpherocylinderSpecies::RunAnalysis() {
  if (sparams_.diffusion_analysis) {
    DiffusionAnalysis();
  }
}

void SpherocylinderSpecies::FinalizeAnalysis() {
  if (sparams_.diffusion_analysis) {
    FinalizeDiffusionAnalysis();
  }
}

void SpherocylinderSpecies::InitDiffusionAnalysis() {
  if (n_members_ == 1) {
    Logger::Warning(
        "Diffusion analysis incompatible with simulations of 1 species member. "
        "Aborting diffusion analysis");
    sparams_.diffusion_analysis = 0;
    return;
  }
  n_samples_ = sparams_.n_diffusion_samples;
  time_ = 0;
  int n_data = params_->n_steps / sparams_.n_posit;
  time_avg_interval_ = n_data / n_samples_;
  if (time_avg_interval_ < 1) {
    Logger::Error("Something went wrong in InitDiffusionAnalysis!");
  }
  pos0_ = new double *[n_members_];
  u0_ = new double *[n_members_];
  for (int i = 0; i < n_members_; ++i) {
    pos0_[i] = new double[params_->n_dim];
    u0_[i] = new double[params_->n_dim];
  }
  vcf_ = new double[time_avg_interval_];
  msd_ = new double[time_avg_interval_];
  vcf_err_ = new double[time_avg_interval_];
  msd_err_ = new double[time_avg_interval_];
  std::fill(msd_, msd_ + time_avg_interval_, 0.0);
  std::fill(vcf_, vcf_ + time_avg_interval_, 0.0);
  std::fill(msd_err_, msd_err_ + time_avg_interval_, 0.0);
  std::fill(vcf_err_, vcf_err_ + time_avg_interval_, 0.0);
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
  for (int i = 0; i < n_members_; ++i) {
    double const *const position0 = members_[i].GetPosition();
    double const *const orientation0 = members_[i].GetOrientation();
    for (int j = 0; j < params_->n_dim; ++j) {
      pos0_[i][j] = position0[j];
      u0_[i][j] = orientation0[j];
    }
  }
}

void SpherocylinderSpecies::CalculateMSD() {
  double avg_sqr_dist = 0.0;
  double avg_sqr_dist_sqr = 0.0;
  for (int i = 0; i < n_members_; ++i) {
    double sqr_diff = 0;
    double const *const position = members_[i].GetPosition();
    for (int j = 0; j < params_->n_dim; ++j) {
      double r_diff = position[j] - pos0_[i][j];
      sqr_diff += SQR(r_diff);
    }
    avg_sqr_dist += sqr_diff;
    avg_sqr_dist_sqr += SQR(sqr_diff);
  }
  avg_sqr_dist /= n_members_;
  avg_sqr_dist_sqr /= n_members_;
  double stdev2 = avg_sqr_dist_sqr - SQR(avg_sqr_dist);
  if (stdev2 < 0) {
    Logger::Error(
        "Something was negative in diffusion analysis when it shouldn't have "
        "been!\n");
  }
  msd_[time_] += avg_sqr_dist / stdev2;
  msd_err_[time_] += 1.0 / stdev2;
}

void SpherocylinderSpecies::CalculateVCF() {
  double avg_udotu0 = 0.0;
  double avg_udotu0_sqr = 0.0;
  for (int i = 0; i < n_members_; ++i) {
    double udotu0 = 0.0;
    double const *const orientation = members_[i].GetOrientation();
    for (int j = 0; j < params_->n_dim; ++j) {
      udotu0 += orientation[j] * u0_[i][j];
    }
    avg_udotu0 += udotu0;
    avg_udotu0_sqr += udotu0 * udotu0;
  }
  avg_udotu0 /= n_members_;
  avg_udotu0_sqr /= n_members_;
  double stdev2 = avg_udotu0_sqr - SQR(avg_udotu0);
  vcf_[time_] += avg_udotu0 / stdev2;
  vcf_err_[time_] += 1.0 / stdev2;
}

void SpherocylinderSpecies::FinalizeDiffusionAnalysis() {
  for (int t = 0; t < time_avg_interval_; ++t) {
    msd_[t] = msd_[t] / msd_err_[t];
    vcf_[t] = vcf_[t] / vcf_err_[t];
    msd_err_[t] = 1.0 / sqrt(n_members_ * msd_err_[t]);
    vcf_err_[t] = 1.0 / sqrt(n_members_ * vcf_err_[t]);
  }
  std::string fname = params_->run_name;
  fname.append("_spherocylinder.diffusion");
  diff_file_.open(fname, std::ios::out);
  diff_file_
      << "length diameter n_dim delta n_steps n_posit n_objs n_samples\n";
  diff_file_ << sparams_.length << " " << sparams_.diameter << " "
             << params_->n_dim << " " << params_->delta << " "
             << params_->n_steps << " " << sparams_.n_posit << " " << n_members_
             << " " << n_samples_ << "\n";
  diff_file_ << "time msd msd_err vcf vcf_err\n";
  diff_file_ << "0.0 0.0 0.0 1.0 0.0\n";
  diff_file_.precision(16);
  diff_file_.setf(std::ios::fixed);
  diff_file_.setf(std::ios::showpoint);
  double midterm = (midstep_ ? 0.5 : 1);
  for (int t = 0; t < time_avg_interval_; ++t) {
    diff_file_ << midterm * (t + 1) * params_->delta * GetNPosit() << " "
               << msd_[t] << " " << msd_err_[t] << " " << vcf_[t] << " "
               << vcf_err_[t] << "\n";
  }
  diff_file_.close();

  for (int i = 0; i < n_members_; ++i) {
    delete pos0_[i];
    delete u0_[i];
  }
  delete pos0_;
  delete u0_;
}
