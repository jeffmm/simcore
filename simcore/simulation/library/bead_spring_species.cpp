#include "bead_spring_species.hpp"
void BeadSpringSpecies::InitAnalysis() {
  time_ = 0;
  // if (params_->bead_spring.diffusion_analysis) {
  // InitDiffusionAnalysis();
  //}
  if (params_->bead_spring.theta_analysis) {
    InitThetaAnalysis();
  }
  if (params_->bead_spring.lp_analysis) {
    InitMse2eAnalysis();
  }
  RunAnalysis();
}

void BeadSpringSpecies::InitMse2eAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_bead_spring.mse2e");
  mse2e_file_.open(fname, std::ios::out);
  mse2e_file_ << "mse2e_analysis_file\n";
  mse2e_file_ << "length diameter bond_length persistence_length driving ndim "
                 "nsteps nspec delta\n";
  auto it = members_.begin();
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  double theory;
  if (params_->n_dim == 2) {
    theory = l * pl * 4.0 - 8.0 * pl * pl * (1 - exp(-0.5 * l / pl));
  } else {
    theory = l * pl * 2.0 - 2.0 * pl * pl * (1 - exp(-l / pl));
  }
  mse2e_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
              << params_->n_dim << " " << params_->n_steps << " " << nspec
              << " " << params_->delta << "\n";
  mse2e_file_ << "num_bead_springs_averaged mse2e_mean mse2e_std_err "
                 "avg_contour_length theory\n";
  mse2e_ = 0;
  mse2e2_ = 0;
  n_samples_ = 0;
  avg_clen_ = 0;
}

void BeadSpringSpecies::InitDiffusionAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_bead_spring.diffusion");
  mse2e_file_.open(fname, std::ios::out);
  mse2e_file_ << "mse2e_analysis_file\n";
  mse2e_file_ << "length diameter bond_length persistence_length driving ndim "
                 "nsteps nspec delta theory\n";
  auto it = members_.begin();
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  double theory;
  if (params_->n_dim == 2) {
    theory = l * pl * 4.0 - 8.0 * pl * pl * (1 - exp(-0.5 * l / pl));
  } else {
    theory = l * pl * 2.0 - 2.0 * pl * pl * (1 - exp(-l / pl));
  }
  mse2e_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
              << params_->n_dim << " " << params_->n_steps << " " << nspec
              << " " << params_->delta << " " << theory << "\n";
  mse2e_file_ << "num_bead_springs_averaged mse2e_mean mse2e_std_err\n";
  mse2e_ = 0.0;
  mse2e2_ = 0.0;
  n_samples_ = 0;
}

void BeadSpringSpecies::InitThetaAnalysis() {
  // TODO Should check to make sure the same lengths, child lengths, persistence
  // lengths, etc are used for each bead_spring in system.
  std::string fname = params_->run_name;
  fname.append("_bead_spring.theta");
  theta_file_.open(fname, std::ios::out);
  theta_file_ << "theta_analysis_file\n";
  theta_file_ << "length diameter bond_length persistence_length "
                 "n_bead_springs n_bonds n_steps n_spec delta n_dim \n";
  double l, cl, pl, dr, d;
  int nbonds;
  int nmembers = members_.size();
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    l = it->GetLength();
    d = it->GetDiameter();
    cl = it->GetBondLength();
    pl = it->GetPersistenceLength();
    dr = it->GetDriving();
    nbonds = it->GetNBonds();
  }
  int nspec = GetNSpec();
  theta_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
              << nmembers << " " << nbonds << " " << params_->n_steps << " "
              << nspec << " " << params_->delta << " " << params_->n_dim << " "
              << "\n";
  theta_file_ << "cos_theta";
  for (int i = 0; i < nbonds - 1; ++i) {
    theta_file_ << " theta_" << i + 1 << i + 2;
  }
  theta_file_ << "\n";
  n_bins_ = 10000;
  int nfil = members_.size();
  theta_histogram_ = new int *[nbonds - 1];
  for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
    theta_histogram_[ibond] = new int[n_bins_];
    for (int ibin = 0; ibin < n_bins_; ++ibin) {
      theta_histogram_[ibond][ibin] = 0;
    }
  }
}

void BeadSpringSpecies::RunAnalysis() {
  // TODO Analyze conformation and ms end-to-end
  if (params_->bead_spring.theta_analysis) {
    if (params_->interaction_flag) {
      std::cout
          << "WARNING! Theta analysis running on interacting bead_springs!\n";
    }
    RunThetaAnalysis();
  }
  if (params_->bead_spring.lp_analysis) {
    RunMse2eAnalysis();
  }
  time_++;
}

void BeadSpringSpecies::RunMse2eAnalysis() {
  // Treat as though we have many spirals for now
  // if ( ! mse2e_file_.is_open()) {
  // early_exit = true;
  // std::cout << " Error! Problem opening file in RunMse2eAnalysis!
  // Exiting.\n";
  //}
  // mse2e_file_ << time_;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    double const *const head_pos = it->GetHeadPosition();
    double const *const tail_pos = it->GetTailPosition();
    double mse2e_temp = 0.0;
    for (int i = 0; i < params_->n_dim; ++i) {
      double temp = (head_pos[i] - tail_pos[i]);
      mse2e_temp += temp * temp;
    }
    mse2e_ += mse2e_temp;
    mse2e2_ += mse2e_temp * mse2e_temp;
    // mse2e_file_ << " " << mse2e ;
    avg_clen_ += it->GetContourLength();
  }
  // mse2e_ /= members_.size();
  // mse2e2_ /= members_.size();
  // mse2e_file_ << "\n";

  n_samples_++;
}

void BeadSpringSpecies::RunThetaAnalysis() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    std::vector<double> const *const thetas = it->GetThetas();
    for (int i = 0; i < (it->GetNBonds() - 1); ++i) {
      int bin_number = (int)floor((1 + (*thetas)[i]) * (n_bins_ / 2));
      if (bin_number == n_bins_) {
        bin_number = n_bins_ - 1;
      } else if (bin_number == -1) {
        bin_number = 0;
      } else if (bin_number > n_bins_ && bin_number < 0) {
        Logger::Error("Something went wrong in RunThetaAnalysis!");
      }
      theta_histogram_[i][bin_number]++;
    }
  }
}

void BeadSpringSpecies::FinalizeAnalysis() {
  if (spiral_file_.is_open()) {
    spiral_file_.close();
  }
  if (theta_file_.is_open()) {
    FinalizeThetaAnalysis();
    theta_file_.close();
  }
  if (mse2e_file_.is_open()) {
    FinalizeMse2eAnalysis();
    mse2e_file_.close();
  }
  if (diffusion_file_.is_open()) {
    // FinalizeDiffusionAnalysis();
    diffusion_file_.close();
  }
}

void BeadSpringSpecies::FinalizeMse2eAnalysis() {
  int num = members_.size();
  mse2e_file_ << num << " ";
  mse2e_ /= n_samples_ * num;
  mse2e2_ /= n_samples_ * num;
  mse2e_file_ << mse2e_ << " ";
  mse2e_file_ << sqrt((mse2e2_ - mse2e_ * mse2e_) / (num * n_samples_)) << " ";
  avg_clen_ /= n_samples_ * num;
  double pl = params_->bead_spring.persistence_length;
  double theory;
  if (params_->n_dim == 2) {
    theory =
        avg_clen_ * pl * 4.0 - 8.0 * pl * pl * (1 - exp(-0.5 * avg_clen_ / pl));
  } else {
    theory = avg_clen_ * pl * 2.0 - 2.0 * pl * pl * (1 - exp(-avg_clen_ / pl));
  }
  mse2e_file_ << avg_clen_ << " " << theory << "\n";
}

void BeadSpringSpecies::FinalizeThetaAnalysis() {
  int nbonds = members_[members_.size() - 1].GetNBonds();
  for (int i = 0; i < n_bins_; ++i) {
    double axis = (2.0 / n_bins_) * i - 1;
    theta_file_ << " " << axis;
    for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
      theta_file_ << " " << theta_histogram_[ibond][i];
    }
    theta_file_ << "\n";
  }

  for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
    delete[] theta_histogram_[ibond];
  }
  delete[] theta_histogram_;
}
