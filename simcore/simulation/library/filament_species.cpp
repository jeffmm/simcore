#include "filament_species.hpp"

FilamentSpecies::FilamentSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::filament);
}
void FilamentSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  fill_volume_ = 0;
  packing_fraction_ = sparams_.packing_fraction;
#ifdef TRACE
  if (packing_fraction_ > 0) {
    Logger::Warning("Simulation run in trace mode with a potentially large "
                    "number of objects in species %s (packing fraction ="
                    " %2.4f)",
                    GetSID()._to_string(), packing_fraction_);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif

  double min_length = 2 * sparams_.min_bond_length;
  if (!sparams_.polydispersity_flag && sparams_.length < min_length) {
    Logger::Warning("Filament length %2.2f is less than minimum filament length"
                    " %2.2f for minimum bond length %2.2f. Setting length to "
                    "minimum length.",
                    sparams_.length, min_length, sparams_.min_bond_length);
    sparams_.length = min_length;
  }
  if (sparams_.perlen_ratio > 0) {
    if (sparams_.polydispersity_flag) {
      Logger::Warning("Persistence length ratio parameter is set with "
                      "polydispersity. Ignoring perlen_ratio and using "
                      "persistence_length parameter instead.");
    } else {
      sparams_.persistence_length = sparams_.length * sparams_.perlen_ratio;
    }
  }
  if (sparams_.spiral_flag && sparams_.spiral_number_fail_condition <= 0) {
    Logger::Warning("Spiral simulation will not end on spiral failure due to"
                    " negative spiral_number_fail_condition");
  }
}

void FilamentSpecies::UpdatePositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  filament_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  filament_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i)
      for (auto it = chunks[i].first; it != chunks[i].second; ++it)
        it->UpdatePosition(midstep_);
  }
#else
  for (filament_iterator it = members_.begin(); it != members_.end(); ++it)
    it->UpdatePosition(midstep_);
#endif

  midstep_ = !midstep_;
}

void FilamentSpecies::Reserve() {
  int max_insert = GetNInsert();
  if (packing_fraction_ > 0) {
    double min_length = 2 * sparams_.min_bond_length;
    double diameter = sparams_.diameter;
    double min_vol = 0;
    if (params_->n_dim == 2) {
      min_vol = diameter * min_length + 0.25 * M_PI * diameter * diameter;
    }
    if (params_->n_dim == 3) {
      min_vol = 0.25 * M_PI * diameter * diameter * min_length +
                1.0 / 6.0 * M_PI * diameter * diameter * diameter;
    }
    max_insert = (int)ceil(packing_fraction_ * space_->volume / min_vol);
  }
  members_.reserve(max_insert);
  Logger::Debug("Reserving memory for %d members in FilamentSpecies",
                max_insert);
}

void FilamentSpecies::AddMember() {
  Species::AddMember();
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ += vol;
    /* if we are still short on volume for the target packing fraction, then
       request more members */
    if (fill_volume_ < packing_fraction_ * space_->volume &&
        members_.size() == sparams_.num) {
      sparams_.num++;
    }
  }
}

void FilamentSpecies::PopMember() {
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ -= vol;
  }
  Species::PopMember();
}

const double FilamentSpecies::GetSpecLength() const {
  if (sparams_.dynamic_instability_flag) {
    return 2 * sparams_.min_bond_length;
  } else {
    return 1.5 * sparams_.min_bond_length;
  }
}
void FilamentSpecies::InitAnalysis() {
  time_ = 0;
  if (sparams_.spiral_flag) {
    InitSpiralAnalysis();
  }
  if (sparams_.theta_analysis) {
    if (params_->interaction_flag) {
      Logger::Warning("Theta analysis running on interacting filaments!");
    }
    InitThetaAnalysis();
  }
  // if (sparams_.crossing_analysis) {
  if (sparams_.crossing_analysis) {
    InitCrossingAnalysis();
  }
  if (sparams_.lp_analysis) {
    InitMse2eAnalysis();
  }
  if (sparams_.global_order_analysis) {
    InitGlobalOrderAnalysis();
  }
  if (params_->polar_order_analysis) {
    InitPolarOrderAnalysis();
  }
  if (sparams_.orientation_corr_analysis) {
    InitOrientationCorrelationAnalysis();
  }
  if (sparams_.flocking_analysis) {
    if (!params_->polar_order_analysis) {
      Logger::Error("Flocking analysis requires polar order analysis to be "
                    "enabled");
    }
    InitFlockingAnalysis();
  }
  RunAnalysis();
  if (params_->in_out_flag) {
    std::string fname = params_->run_name;
    fname.append("_filament.in_out");
    in_out_file_.open(fname, std::ios::out);
    in_out_file_ << "perlen: " << sparams_.perlen_ratio << "\n";
    in_out_file_ << "umax: " << params_->soft_potential_mag << "\n";
    in_out_file_ << "angle_in: ";
    if (members_.size() < 2) {
      Logger::Error(
          "Error in in-out analysis in FilamentSpecies::InitAnalysis");
    }
    double u1[3] = {0, 0, 0};
    double u2[3] = {0, 0, 0};
    members_[0].GetAvgOrientation(u1);
    members_[1].GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim, u1, u2);
    in_out_file_ << acos(dp) << "\n";
  }
}

void FilamentSpecies::InitGlobalOrderAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.global_order");
  global_order_file_.open(fname, std::ios::out);
  global_order_file_ << "global_order_analysis_file\n";
  global_order_file_
      << "time polar_order_x polar_order_y polar_order_z nematic_order_xx "
         "nematic_order_xy nematic_order_xz nematic_order_yx nematic_order_yy "
         "nematic_order_yz nematic_order_zx nematic_order_zy nematic_order_zz "
         "spiral_order signed_spiral_order n_spooling "
         "avg_spool_spiral_number\n";
  nematic_order_tensor_ = new double[9];
  polar_order_vector_ = new double[3];
  std::fill(nematic_order_tensor_, nematic_order_tensor_ + 9, 0.0);
  std::fill(polar_order_vector_, polar_order_vector_ + 3, 0.0);
}

void FilamentSpecies::InitPolarOrderAnalysis() {
  n_bins_1d_ = params_->polar_order_n_bins;
  // Ensure n_bins_1d_ is even, to avoid headaches
  if (n_bins_1d_ % 2 != 0) {
    n_bins_1d_++;
  }
  n_bins_ = SQR(n_bins_1d_);
  polar_order_histogram_ = new int[n_bins_];
  std::fill(polar_order_histogram_, polar_order_histogram_ + n_bins_, 0.0);
  contact_cut_ = params_->polar_order_contact_cutoff;
  contact_bin_width_ = contact_cut_ / n_bins_1d_;
  polar_bin_width_ = 2.0 / n_bins_1d_;
  std::string fname = params_->run_name;
  fname.append("_filament.polar_order");
  polar_order_file_.open(fname, std::ios::out);
  fname = params_->run_name;
  fname.append("_filament.polar_order_avg");
  polar_order_avg_file_.open(fname, std::ios::out);
  polar_order_avg_file_ << "polar_order_avg_file\n";
  polar_order_file_ << "contact_number local_polar_order\n";
  polar_order_avg_file_ << "time avg_polar_order avg_contact_number\n";
}

void FilamentSpecies::RunPolarOrderAnalysis() {
  std::vector<double> po;
  std::vector<double> cn;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->GetPolarOrders(&po);
    it->GetContactNumbers(&cn);
  }
  if (po.size() != cn.size()) {
    Logger::Error(
        "Number of polar order parameters and contact numbers not equal");
  }
  po_avg_ = 0;
  cn_avg_ = 0;
  for (int i = 0; i < po.size(); ++i) {
    po_avg_ += po[i];
    cn_avg_ += cn[i];
    if (cn[i] > contact_cut_) {
      continue;
    }
    // polar_order_file_ << cn[i] << " " << po[i] <<"\n";
    int x = (int)(floor(cn[i] / contact_bin_width_));
    int y = (int)(floor((po[i] + 1) / polar_bin_width_));
    if (y == n_bins_1d_)
      y = n_bins_1d_ - 1;
    if (x == n_bins_1d_)
      x = n_bins_1d_ - 1;
    if (y == -1)
      y = 0;
    if (x == -1)
      x = 0;
    if (y < 0 || x < 0 || y > n_bins_1d_ - 1 || x > n_bins_1d_ - 1) {
      std::cout << cn[i] << " " << po[i] << "\n";
      Logger::Error("Out of range in RunPolarOrderAnalysis");
    }
    polar_order_histogram_[n_bins_1d_ * y + x]++;
  }
  po_avg_ /= po.size();
  cn_avg_ /= cn.size();
  polar_order_avg_file_ << time_ << " " << po_avg_ << " " << cn_avg_ << "\n";
}

void FilamentSpecies::InitOrientationCorrelationAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.orientation_corr");
  orientation_corr_file_.open(fname, std::ios::out);
  orientation_corr_n_steps_ = sparams_.orientation_corr_n_steps;
  orientation_corr_file_ << "orientation_corr_analysis_file, n_filaments = "
                         << n_members_
                         << ", n_bonds = " << members_[0].GetNBonds()
                         << ", n_avg_steps = " << orientation_corr_n_steps_
                         << "\n";
  orientation_corr_file_ << "time orientation_corr_avg orientation_corr_sem\n";
}

void FilamentSpecies::RunOrientationCorrelationAnalysis() {
  if (time_ % orientation_corr_n_steps_ == 0) {
    for (auto it = members_.begin(); it != members_.end(); ++it) {
      it->ZeroOrientationCorrelations();
    }
    orientation_corr_file_ << "0 1 0\n";
    return;
  }
  double avg = 0;
  double sem = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    std::pair<double, double> avg_var = it->GetAvgOrientationCorrelation();
    avg += avg_var.first;
    sem += avg_var.second;
  }
  avg /= n_members_;
  sem = sem / n_members_; // This calculates the pooled variance
  sem =
      sqrt(sem / (n_members_ * members_[0].GetNBonds())); // """ standard error
  orientation_corr_file_ << time_ % orientation_corr_n_steps_ << " " << avg
                         << " " << sem << "\n";
}

void FilamentSpecies::FinalizeOrientationCorrelationAnalysis() {}

void FilamentSpecies::InitCrossingAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.crossing");
  crossing_file_.open(fname, std::ios::out);
  crossing_file_ << "crossing_analysis\n";
  crossing_file_ << "sp lp dr\n";
  crossing_file_ << params_->soft_potential_mag << " " << sparams_.perlen_ratio
                 << " " << sparams_.driving_factor << "\n";
  double avg_u[3] = {0, 0, 0};
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->GetAvgOrientation(avg_u);
    crossing_file_ << acos(avg_u[1]) << " ";
  }
  crossing_file_ << "\n";
}

void FilamentSpecies::RunCrossingAnalysis() {}

void FilamentSpecies::FinalizeCrossingAnalysis() {
  double avg_pos[3] = {0, 0, 0};
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->GetAvgPosition(avg_pos);
    crossing_file_ << (avg_pos[0] > 2 ? 1 : 0) << " ";
  }
  crossing_file_ << "\n";
}

void FilamentSpecies::InitMse2eAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.mse2e");
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
  // if (params_->n_dim == 2) {
  // theory = l * pl * 4.0 - 8.0 * pl * pl * (1 - exp(-0.5 * l / pl));
  //} else {
  theory = l * pl * 2.0 - 2.0 * pl * pl * (1 - exp(-l / pl));
  //}
  mse2e_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
              << params_->n_dim << " " << params_->n_steps << " " << nspec
              << " " << params_->delta << " " << theory << "\n";
  mse2e_file_ << "num_filaments_averaged mse2e_mean mse2e_std_err\n";
  mse2e_ = 0.0;
  mse2e2_ = 0.0;
  n_samples_ = 0;
}

void FilamentSpecies::InitSpiralAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.spiral");
  spiral_file_.open(fname, std::ios::out);
  spiral_file_ << "spiral_analysis_file\n";
  spiral_file_ << "length diameter bond_length persistence_length driving "
                  "nsteps nspec delta\n";
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    double l = it->GetLength();
    double d = it->GetDiameter();
    double cl = it->GetBondLength();
    double pl = it->GetPersistenceLength();
    double dr = it->GetDriving();
    double nspec = GetNSpec();
    spiral_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
                 << params_->n_steps << " " << nspec << " " << params_->delta
                 << "\n";
  }
  spiral_file_ << "time angle_sum E_bend tip_z_proj spiral_number head_pos_x "
                  "head_pos_y tail_pos_x tail_pos_y\n";
}

void FilamentSpecies::InitThetaAnalysis() {
  // TODO Should check to make sure the same lengths, child lengths, persistence
  // lengths, etc are used for each filament in system.
  std::string fname = params_->run_name;
  fname.append("_filament.theta");
  theta_file_.open(fname, std::ios::out);
  theta_file_ << "theta_analysis_file\n";
  theta_file_
      << "length diameter bond_length persistence_length driving n_filaments "
         "n_bonds n_steps n_spec delta n_dim \n";
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
              << nspec << " " << params_->delta << " " << params_->n_dim
              << "\n";
  theta_file_ << "cos_theta";
  for (int i = 0; i < nbonds - 1; ++i) {
    theta_file_ << " theta_" << i + 1 << i + 2;
  }
  theta_file_ << "\n";
  n_bins_ = 10000;
  theta_histogram_ = new int *[nbonds - 1];
  for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
    theta_histogram_[ibond] = new int[n_bins_];
    for (int ibin = 0; ibin < n_bins_; ++ibin) {
      theta_histogram_[ibond][ibin] = 0;
    }
  }
}

void FilamentSpecies::RunAnalysis() {
  if (sparams_.spiral_flag) {
    RunSpiralAnalysis();
  }
  // TODO Analyze conformation and ms end-to-end
  if (sparams_.theta_analysis) {
    RunThetaAnalysis();
  }
  if (sparams_.lp_analysis) {
    RunMse2eAnalysis();
  }
  if (sparams_.global_order_analysis) {
    RunGlobalOrderAnalysis();
  }
  if (params_->polar_order_analysis) {
    RunPolarOrderAnalysis();
  }
  if (sparams_.orientation_corr_analysis) {
    RunOrientationCorrelationAnalysis();
  }
  if (sparams_.flocking_analysis) {
    RunFlockingAnalysis();
  }
  time_++;
}

void FilamentSpecies::RunGlobalOrderAnalysis() {
  double sn_tot = 0.0;
  double sn_mag = 0.0;
  double sn_spools = 0;
  int n_spooling = 0;
  double sn;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->CalculateSpiralNumber();
    sn = it->GetSpiralNumber();
    it->GetPolarOrder(polar_order_vector_);
    it->GetNematicOrder(nematic_order_tensor_);
    sn_mag += ABS(sn);
    sn_tot += sn;
    if (sn > 0.7) {
      sn_spools += sn;
      n_spooling++;
    }
  }
  sn_mag /= n_members_;
  sn_tot /= n_members_;
  if (n_spooling > 0) {
    sn_spools /= n_spooling;
  } else {
    sn_spools = 0;
  }
  for (int i = 0; i < 3; ++i) {
    polar_order_vector_[i] /= n_members_;
  }
  for (int i = 0; i < 9; ++i) {
    nematic_order_tensor_[i] /= n_members_;
  }
  if (global_order_file_.is_open()) {
    global_order_file_ << time_ << " " << polar_order_vector_[0] << " "
                       << polar_order_vector_[1] << " "
                       << polar_order_vector_[2] << " "
                       << nematic_order_tensor_[0] << " "
                       << nematic_order_tensor_[1] << " "
                       << nematic_order_tensor_[2] << " "
                       << nematic_order_tensor_[3] << " "
                       << nematic_order_tensor_[4] << " "
                       << nematic_order_tensor_[5] << " "
                       << nematic_order_tensor_[6] << " "
                       << nematic_order_tensor_[7] << " "
                       << nematic_order_tensor_[8] << " " << sn_mag << " "
                       << sn_tot << " " << n_spooling << " " << sn_spools
                       << "\n";
  } else {
    Logger::Error("Problem opening file in RunGlobalOrderAnalysis!");
  }
}

void FilamentSpecies::RunSpiralAnalysis() {
  // Treat as though we have many spirals for now
  double tip_z;
  auto it = members_.begin();
  e_bend_ = tot_angle_ = 0;
  double length = it->GetLength();
  double plength = it->GetPersistenceLength();
  double clength = it->GetBondLength();
  double e_zero = length * plength / (clength * clength);
  it->CalculateSpiralNumber();
  double spiral_number = it->GetSpiralNumber();
  std::vector<double> const *const thetas = it->GetThetas();
  for (int i = 0; i < thetas->size(); ++i) {
    tot_angle_ += acos((*thetas)[i]);
    e_bend_ += (*thetas)[i];
  }
  // record energy relative to the bending "zero energy" (straight rod)
  e_bend_ = e_zero - e_bend_ * plength / clength;
  tip_z = it->GetTipZ();
  double const *const head_pos = it->GetHeadPosition();
  double const *const tail_pos = it->GetTailPosition();
  if (spiral_file_.is_open()) {
    spiral_file_ << time_ << " " << tot_angle_ << " " << e_bend_ << " " << tip_z
                 << " " << spiral_number << " " << head_pos[0] << " "
                 << head_pos[1] << " " << tail_pos[0] << " " << tail_pos[1]
                 << "\n";
  } else {
    Logger::Error("Problem opening file in RunSpiralAnalysis!");
  }
}

void FilamentSpecies::RunMse2eAnalysis() {
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
  }
  // mse2e_ /= members_.size();
  // mse2e2_ /= members_.size();
  // mse2e_file_ << "\n";
  n_samples_++;
}

void FilamentSpecies::RunThetaAnalysis() {
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

void FilamentSpecies::FinalizeAnalysis() {
  if (spiral_file_.is_open()) {
    spiral_file_.close();
  }
  if (theta_file_.is_open()) {
    FinalizeThetaAnalysis();
    theta_file_.close();
  }
  if (crossing_file_.is_open()) {
    FinalizeCrossingAnalysis();
    crossing_file_.close();
  }
  if (mse2e_file_.is_open()) {
    FinalizeMse2eAnalysis();
    mse2e_file_.close();
  }
  if (global_order_file_.is_open()) {
    FinalizeGlobalOrderAnalysis();
    global_order_file_.close();
  }
  if (polar_order_file_.is_open()) {
    FinalizePolarOrderAnalysis();
    polar_order_file_.close();
    polar_order_avg_file_.close();
  }
  if (orientation_corr_file_.is_open()) {
    FinalizeOrientationCorrelationAnalysis();
    orientation_corr_file_.close();
  }
  if (flock_file_.is_open()) {
    FinalizeFlockingAnalysis();
    flock_file_.close();
  }
  if (in_out_file_.is_open()) {
    in_out_file_ << "angle_out: ";
    if (members_.size() < 2) {
      Logger::Error(
          "Error in in-out analysis in FilamentSpecies::FinalizeAnalysis");
    }
    double u1[3] = {0, 0, 0};
    double u2[3] = {0, 0, 0};
    members_[0].GetAvgOrientation(u1);
    members_[1].GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim, u1, u2);
    in_out_file_ << acos(dp) << "\n";
    in_out_file_.close();
  }
}

void FilamentSpecies::FinalizeGlobalOrderAnalysis() {}

void FilamentSpecies::FinalizePolarOrderAnalysis() {
  /* In order to avoid overcounting cases where there were no local interactors
   * to count for local polar order, I am going to smooth the bin representing
   * (0,0) in the histogram, by averaging vertically along the y-axis */
  int avg_bin = (int)floor(
      0.5 * (polar_order_histogram_[(n_bins_1d_ / 2 - 1) * n_bins_1d_] +
             polar_order_histogram_[(n_bins_1d_ / 2 + 1) * n_bins_1d_]));
  polar_order_histogram_[(n_bins_1d_ / 2) * n_bins_1d_] = avg_bin;
  for (int i = 0; i < n_bins_1d_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      polar_order_file_
          << polar_order_histogram_[(n_bins_1d_ - 1 - i) * n_bins_1d_ + j]
          << " ";
    }
    polar_order_file_ << "\n";
  }
}

// void FilamentSpecies::FinalizeLocalOrderAnalysis() {
// if (params_->local_structure_average) {
// WriteLocalOrderData();
//}
//}

void FilamentSpecies::FinalizeMse2eAnalysis() {
  int num = members_.size();
  mse2e_file_ << num << " ";
  mse2e_ /= n_samples_ * num;
  mse2e2_ /= n_samples_ * num;
  mse2e_file_ << mse2e_ << " ";
  mse2e_file_ << sqrt((mse2e2_ - mse2e_ * mse2e_) / (num * n_samples_)) << "\n";
}

void FilamentSpecies::FinalizeThetaAnalysis() {
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

void FilamentSpecies::InitFlockingAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.flock");
  flock_file_.open(fname, std::ios::out);
  flock_file_ << "flock_analysis_file\n";
  flock_file_ << "length diameter bond_length persistence_length umax driving "
                 "nsteps nspec delta\n";
  auto it = members_.begin();
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  flock_file_ << l << " " << d << " " << cl << " " << pl << " "
              << params_->soft_potential_mag << " " << dr << " "
              << params_->n_steps << " " << nspec << " " << params_->delta
              << "\n";
  flock_file_ << "time n_flocking n_exterior n_interior n_joined n_left";
  for (int i = 0; i < members_.size(); ++i) {
    flock_file_ << " fil" << i;
  }
  flock_file_ << "\n";
  flock_states_ = new int[members_.size()];
}

void FilamentSpecies::RunFlockingAnalysis() {
  int n_flocking = 0;
  int n_interior = 0;
  int n_exterior = 0;
  int n_joined = 0;
  int n_left = 0;
  for (int i = 0; i < members_.size(); ++i) {

    members_[i].CheckFlocking();
    int flock_type = members_[i].GetFlockType();
    int flock_change_state = members_[i].GetFlockChangeState();
    flock_states_[i] = flock_type;
    if (flock_type) {
      n_flocking++;
      if (flock_type == 1) {
        n_exterior++;
      } else if (flock_type == 2) {
        n_interior++;
      } else {
        Logger::Warning("Flock type not recognized in "
                        "FilamentSpecies::RunFlockingAnalysis");
      }
    }
    if (flock_change_state) {
      if (flock_change_state == 1) {
        n_joined++;
      } else if (flock_change_state == 2) {
        n_left++;
      } else {
        Logger::Warning("Flock change state not recognized in "
                        "FilamentSpecies::RunFlockingAnalysis");
      }
    }
  }
  if (flock_file_.is_open()) {
    flock_file_ << time_ << " " << n_flocking << " " << n_exterior << " "
                << n_interior << " " << n_joined << " " << n_left;
    for (int i=0; i<members_.size(); ++i) {
      flock_file_ << " " << flock_states_[i];
    }
    flock_file_ << "\n";
  } else {
    Logger::Error("Problem opening file in RunFlockingAnalysis!");
  }
}

void FilamentSpecies::FinalizeFlockingAnalysis() { delete[] flock_states_; }
