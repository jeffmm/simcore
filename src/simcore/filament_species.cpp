#include "filament.h"

void FilamentSpecies::InitAnalysis() {
  time_ = 0;
  if (params_->filament.spiral_flag) {
    InitSpiralAnalysis();
  }
  if (params_->filament.theta_analysis) {
    if (params_->interaction_flag) {
      std::cout << "WARNING! Theta analysis running on interacting filaments!\n";
    }
    InitThetaAnalysis();
  }
  //if (params_->filament.crossing_analysis) {
  if (params_->filament.crossing_analysis) {
    InitCrossingAnalysis();
  }
  if (params_->filament.lp_analysis) {
    InitMse2eAnalysis();
  }
  if (params_->filament.global_order_analysis) {
    InitGlobalOrderAnalysis();
  }
  if (params_->polar_order_analysis) {
    InitPolarOrderAnalysis();
  }
  if (params_->filament.orientation_corr_analysis) {
    InitOrientationCorrelationAnalysis();
  }
  if (params_->filament.flocking_analysis) {
    if (!params_->polar_order_analysis) {
      warning("Flocking analysis enabled polar order analysis in FilamentSpecies::InitAnalysis\n");
      params_->polar_order_analysis = 1;
    }
    InitFlockingAnalysis();
  }
  RunAnalysis();
  if (params_->in_out_flag) {
    std::string fname = params_->run_name;
    fname.append("_filament.in_out");
    in_out_file_.open(fname, std::ios::out);
    in_out_file_ << "perlen: " << params_->filament.perlen_ratio << "\n";
    in_out_file_ << "umax: " << params_->soft_potential_mag << "\n";
    in_out_file_ << "angle_in: ";
    if (members_.size() < 2) {
      error_exit("Error in in-out analysis in FilamentSpecies::InitAnalysis");
    }
    double u1[3] = {0,0,0};
    double u2[3] = {0,0,0};
    members_[0].GetAvgOrientation(u1);
    members_[1].GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim,u1,u2);
    in_out_file_ << acos(dp) <<"\n";
  }
}

void FilamentSpecies::InitGlobalOrderAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.global_order");
  global_order_file_.open(fname, std::ios::out);
  global_order_file_ << "global_order_analysis_file\n";
  global_order_file_ << "time polar_order_x polar_order_y polar_order_z nematic_order_xx nematic_order_xy nematic_order_xz nematic_order_yx nematic_order_yy nematic_order_yz nematic_order_zx nematic_order_zy nematic_order_zz spiral_order signed_spiral_order n_spooling avg_spool_spiral_number\n";
  nematic_order_tensor_ = new double[9];
  polar_order_vector_ = new double[3];
  std::fill(nematic_order_tensor_, nematic_order_tensor_+9, 0.0);
  std::fill(polar_order_vector_, polar_order_vector_+3, 0.0);
}



void FilamentSpecies::InitPolarOrderAnalysis() {
  n_bins_1d_ = params_->polar_order_n_bins;
  // Ensure n_bins_1d_ is even, to avoid headaches
  if (n_bins_1d_%2 != 0) {
    n_bins_1d_++;
  }
  n_bins_ = SQR(n_bins_1d_);
  polar_order_histogram_ = new int[n_bins_];
  std::fill(polar_order_histogram_,polar_order_histogram_+n_bins_,0.0);
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
  for (auto it=members_.begin(); it!= members_.end(); ++it) {
    it->GetPolarOrders(&po);
    it->GetContactNumbers(&cn);
  }
  if (po.size() != cn.size()) {
    error_exit("Number of polar order parameters and contact numbers not equal");
  }
  po_avg_ = 0;
  cn_avg_ = 0;
  for (int i=0; i<po.size(); ++i) {
    po_avg_ += po[i];
    cn_avg_ += cn[i];
    if (cn[i] > contact_cut_) {
      continue;
    }
    //polar_order_file_ << cn[i] << " " << po[i] <<"\n";
    int x = (int) (floor(cn[i]/contact_bin_width_));
    int y = (int) (floor((po[i]+1)/polar_bin_width_));
    if (y == n_bins_1d_) y = n_bins_1d_-1;
    if (x == n_bins_1d_) x = n_bins_1d_-1;
    if (y == -1) y = 0;
    if (x == -1) x = 0;
    if (y<0 || x<0 || y>n_bins_1d_-1 || x>n_bins_1d_-1) {
      std::cout << cn[i] << " " << po[i] << "\n";
      error_exit("Out of range in RunPolarOrderAnalysis");
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
  orientation_corr_n_steps_ = params_->filament.orientation_corr_n_steps;
  orientation_corr_file_ << "orientation_corr_analysis_file, n_filaments = " << n_members_ << ", n_bonds = " << members_[0].GetNBonds() << ", n_avg_steps = " << orientation_corr_n_steps_ << "\n";
  orientation_corr_file_ << "time orientation_corr_avg orientation_corr_sem\n";
}

void FilamentSpecies::RunOrientationCorrelationAnalysis() {
  if (time_%orientation_corr_n_steps_ == 0) {
    for (auto it=members_.begin(); it!=members_.end(); ++it) {
      it->ZeroOrientationCorrelations();
    }
    orientation_corr_file_ << "0 1 0\n";
    return;
  }
  double avg = 0;
  double sem = 0;
  for (auto it=members_.begin(); it!= members_.end(); ++it) {
    std::pair<double, double> avg_var = it->GetAvgOrientationCorrelation();
    avg += avg_var.first;
    sem += avg_var.second;
  }
  avg/=n_members_;
  sem = sem / n_members_; // This calculates the pooled variance
  sem = sqrt(sem/(n_members_*members_[0].GetNBonds())); // """ standard error
  orientation_corr_file_ << time_%orientation_corr_n_steps_ << " " << avg << " " << sem << "\n";
}

void FilamentSpecies::FinalizeOrientationCorrelationAnalysis() {

}

void FilamentSpecies::InitCrossingAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.crossing");
  crossing_file_.open(fname, std::ios::out);
  crossing_file_ << "crossing_analysis\n";
  crossing_file_ << "sp lp dr\n";
  crossing_file_ << params_->soft_potential_mag << " " << params_->filament.perlen_ratio << " " << params_->filament.driving_factor <<"\n";
  double avg_u[3] = {0,0,0};
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    it->GetAvgOrientation(avg_u);
    crossing_file_ << acos(avg_u[1]) << " ";
  }
  crossing_file_ << "\n";
}

void FilamentSpecies::RunCrossingAnalysis() {

}

void FilamentSpecies::FinalizeCrossingAnalysis() {
  double avg_pos[3] = {0,0,0};
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
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
  mse2e_file_ << "length diameter bond_length persistence_length driving ndim nsteps nspec delta theory\n";
  auto it=members_.begin();
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  double theory;
  if (params_->n_dim == 2) {
    theory = l * pl * 4.0 - 8.0 * pl * pl * (1-exp(-0.5*l/pl));
  }
  else {
    theory = l * pl * 2.0 - 2.0 * pl * pl * (1-exp(-l/pl));
  }
  mse2e_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " " << params_->n_dim << " " << params_->n_steps << " " << nspec << " " << params_->delta << " " << theory << "\n";
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
  spiral_file_ << "length diameter bond_length persistence_length driving nsteps nspec delta\n";
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    double l = it->GetLength();
    double d = it->GetDiameter();
    double cl = it->GetBondLength();
    double pl = it->GetPersistenceLength();
    double dr = it->GetDriving();
    double nspec = GetNSpec();
    spiral_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " " << params_->n_steps << " " << nspec << " " << params_->delta << "\n";
  }
  spiral_file_ << "time angle_sum E_bend tip_z_proj spiral_number head_pos_x head_pos_y tail_pos_x tail_pos_y\n";
}

void FilamentSpecies::InitThetaAnalysis() {
  // TODO Should check to make sure the same lengths, child lengths, persistence lengths, etc are used for each filament in system.
  std::string fname = params_->run_name;
  fname.append("_filament.theta");
  theta_file_.open(fname, std::ios::out);
  theta_file_ << "theta_analysis_file\n";
  theta_file_ << "length diameter bond_length persistence_length driving n_filaments n_bonds n_steps n_spec delta n_dim metric_forces\n";
  double l, cl, pl, dr, d;
  int nbonds;
  int nmembers = members_.size();
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    l = it->GetLength();
    d = it->GetDiameter();
    cl = it->GetBondLength();
    pl = it->GetPersistenceLength();
    dr = it->GetDriving();
    nbonds = it->GetNBonds();
  }
  int nspec = GetNSpec();
  theta_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " " << nmembers << " " << nbonds << " " << params_->n_steps << " " << nspec << " " << params_->delta << " " << params_->n_dim << " " << params_->filament.metric_forces << "\n";
  theta_file_ << "cos_theta";
  for (int i=0; i<nbonds-1; ++i) {
    theta_file_ << " theta_" << i+1 << i+2;
  }
  theta_file_ << "\n";
  n_bins_ = 10000;
  int nfil = members_.size();
  theta_histogram_ = new int *[nbonds-1];
  for (int ibond=0; ibond<nbonds-1; ++ibond) {
    theta_histogram_[ibond] = new int[n_bins_];
    for (int ibin=0;ibin<n_bins_;++ibin) {
      theta_histogram_[ibond][ibin] = 0;
    }
  }
}

void FilamentSpecies::RunAnalysis() {
  if (params_->filament.spiral_flag) {
    RunSpiralAnalysis();
  }
  // TODO Analyze conformation and ms end-to-end
  if (params_->filament.theta_analysis) {
    RunThetaAnalysis();
  }
  if (params_->filament.lp_analysis) {
    RunMse2eAnalysis();
  }
  if (params_->filament.global_order_analysis) {
    RunGlobalOrderAnalysis();
  }
  if (params_->polar_order_analysis) {
    RunPolarOrderAnalysis();
  }
  if (params_->filament.orientation_corr_analysis) {
    RunOrientationCorrelationAnalysis();
  }
  if (params_->filament.flocking_analysis) {
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
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
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
  }
  else {
    sn_spools = 0;
  }
  for (int i=0; i<3; ++i) {
    polar_order_vector_[i]/=n_members_;
  }
  for (int i=0; i<9; ++i) {
    nematic_order_tensor_[i]/=n_members_;
  }
  if (global_order_file_.is_open()) {
    global_order_file_ << time_ << " " << polar_order_vector_[0] << " " << polar_order_vector_[1] << " " << polar_order_vector_[2] << " " << nematic_order_tensor_[0] << " " << nematic_order_tensor_[1] << " " << nematic_order_tensor_[2] << " " << nematic_order_tensor_[3] << " " << nematic_order_tensor_[4] << " " << nematic_order_tensor_[5] << " " << nematic_order_tensor_[6] << " " << nematic_order_tensor_[7] << " " << nematic_order_tensor_[8] << " " << sn_mag << " " << sn_tot << " " << n_spooling << " " << sn_spools << "\n";
  }
  else {
    early_exit = true;
    std::cout << "ERROR: Problem opening file in RunGlobalOrderAnalysis! Exiting.\n";
  }
}

//void FilamentSpecies::RunLocalOrderAnalysis() {
  //// Calculate extra-filament local orientation correlations
  //// All-pairs minimum distance
  //// Want a local density within a radius of 1.5*length
  //double rcut2 = SQR(0.5*n_bins_1d_*params_->local_structure_bin_width);
  ////double rcut2 = SQR(1.5)*params_->filament.length*params_->filament.length;
  //for (auto it=members_.begin(); it!= members_.end(); ++it) {
    //std::vector<Interaction*> * ixs = it->GetInteractions();
    //double const * const r1 = it->GetPosition();
    //double const * const u1 = it->GetOrientation();
    //// Store the MeshIDs of all objects within range of filament_i
    //std::set<int> ix_mids;
    //for (auto ix = ixs->begin(); ix!=ixs->end(); ++ix) {
      //if ((*ix)->mids.first != it->GetMeshID()) {
        //ix_mids.insert((*ix)->mids.first);
      //}
      //else if ((*ix)->mids.second != it->GetMeshID()) {
        //ix_mids.insert((*ix)->mids.second);
      //}
    //}
    //for (auto jt=members_.begin(); jt!=members_.end(); ++jt) {
      //// Avoid intra-filament correlations for now
      //if (it->GetMeshID() == jt->GetMeshID()) continue;
      //// Check if filament_j is within range
      //if (ix_mids.find(jt->GetMeshID()) == ix_mids.end()) continue;

      //std::cout << " Calculating interaction\n";
      //double const * const r2 = jt->GetPosition();
      //double r_diff[3] = {0,0,0};
      //double r_diff_mag = 0;
      //for (int i=0; i<params_->n_dim; ++i) {
        //r_diff[i] = r2[i] - r1[i];
        //r_diff_mag += r_diff[i]*r_diff[i];
      //}
      ////if (r_diff_mag > 4*params_->filament.length*params_->filament.length) continue;
      ////mindist_.ObjectObject(&(*it),&(*jt),&ix);
      ////if (ix.dr_mag2 > rcut2) continue;
      ////double const * const u2 = jt->GetOrientation();
      ////double dp = dot_product(params_->n_dim,u1,u2);
      ////double temp[3];
      ////cross_product(u1,u2,temp,3);
      ////double relative_angle = SIGNOF(temp[2]) * acos(dp);

      //// Here, we want to use Bresenham's line drawing algorithm to populate the histograms

    //}
  //}
  //// If we are doing a time average, wait to write local order data.
  //if (params_->local_structure_average == 0 ) {
  //// Write local order data
    //WriteLocalOrderData();
  //}
//}

//void FilamentSpecies::WriteLocalOrderData() {
  //if (local_order_file_.is_open()) {
    //local_order_file_ << time_ << " " << n_bins_ << " ";
    //for (int i=0;i<n_bins_1d_;++i) {
      //for (int j=0; j<n_bins_1d_; ++j) {
        //local_order_file_ << pdf_histogram_[i][j];
      //}
    //}
    //local_order_file_ << "\n";
  //}
  //else {
    //early_exit = true;
    //std::cout << "ERROR: Problem opening file in RunLocalOrderAnalysis! Exiting.\n";
  //}

//}


void FilamentSpecies::RunSpiralAnalysis() {
  // Treat as though we have many spirals for now
  double tip_z;
  auto it=members_.begin();
  e_bend_ = tot_angle_ = 0;
  double length = it->GetLength();
  double plength = it->GetPersistenceLength();
  double clength = it->GetBondLength();
  double e_zero = length * plength / (clength * clength);
  it->CalculateSpiralNumber();
  double spiral_number = it->GetSpiralNumber();
  std::vector<double> const * const thetas = it->GetThetas();
  for (int i=0; i<thetas->size(); ++i) {
    tot_angle_ += acos((*thetas)[i]);
    e_bend_ += (*thetas)[i];
  }
  // record energy relative to the bending "zero energy" (straight rod)
  e_bend_ = e_zero - e_bend_ * plength / clength;
  tip_z = it->GetTipZ();
  double const * const head_pos = it->GetHeadPosition();
  double const * const tail_pos = it->GetTailPosition();
  if (spiral_file_.is_open()) {
    spiral_file_ << time_ << " " << tot_angle_ << " " << e_bend_ << " " << tip_z << " " << spiral_number << " " << head_pos[0] << " " << head_pos[1] << " " << tail_pos[0] << " " << tail_pos[1] << "\n";
  }
  else {
    early_exit = true;
    std::cout << "ERROR: Problem opening file in RunSpiralAnalysis! Exiting.\n";
  }
}

void FilamentSpecies::RunMse2eAnalysis() {
  // Treat as though we have many spirals for now
  //if ( ! mse2e_file_.is_open()) {
    //early_exit = true;
    //std::cout << " Error! Problem opening file in RunMse2eAnalysis! Exiting.\n";
  //}
  //mse2e_file_ << time_;
  for (auto it=members_.begin(); it!= members_.end(); ++it) {
    double const * const head_pos = it->GetHeadPosition();
    double const * const tail_pos = it->GetTailPosition();
    double mse2e_temp = 0.0;
    for (int i=0; i<params_->n_dim; ++i) {
      double temp = (head_pos[i] - tail_pos[i]);
      mse2e_temp += temp*temp;
    }
    mse2e_ += mse2e_temp;
    mse2e2_ += mse2e_temp*mse2e_temp;
    //mse2e_file_ << " " << mse2e ;
  }
  //mse2e_ /= members_.size();
  //mse2e2_ /= members_.size();
  //mse2e_file_ << "\n";
  n_samples_++;
}

void FilamentSpecies::RunThetaAnalysis() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    std::vector<double> const * const thetas = it->GetThetas();
    for (int i=0; i<(it->GetNBonds()-1); ++i) {
      int bin_number = (int) floor( (1 + (*thetas)[i]) * (n_bins_/2) );
      if (bin_number == n_bins_) {
        bin_number = n_bins_-1;
      }
      else if (bin_number == -1) {
        bin_number = 0;
      }
      else if (bin_number > n_bins_ && bin_number < 0) {
        error_exit("Something went wrong in RunThetaAnalysis!");
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
      error_exit("Error in in-out analysis in FilamentSpecies::FinalizeAnalysis");
    }
    double u1[3] = {0,0,0};
    double u2[3] = {0,0,0};
    members_[0].GetAvgOrientation(u1);
    members_[1].GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim,u1,u2);
    in_out_file_ << acos(dp) <<"\n";
    in_out_file_.close();
  }
}

void FilamentSpecies::FinalizeGlobalOrderAnalysis() {
}

void FilamentSpecies::FinalizePolarOrderAnalysis() {
  /* In order to avoid overcounting cases where there were no local interactors to
   * count for local polar order, I am going to smooth the bin representing (0,0)
   * in the histogram, by averaging vertically along the y-axis */
  int avg_bin = (int) floor(0.5*(polar_order_histogram_[(n_bins_1d_/2-1)*n_bins_1d_]
        + polar_order_histogram_[(n_bins_1d_/2+1)*n_bins_1d_]));
  polar_order_histogram_[(n_bins_1d_/2)*n_bins_1d_] = avg_bin;
  for (int i=0;i<n_bins_1d_; ++i) {
    for (int j=0;j<n_bins_1d_;++j) {
      polar_order_file_ << polar_order_histogram_[(n_bins_1d_-1-i)*n_bins_1d_+j] << " ";
    }
    polar_order_file_ << "\n";
  }
}

//void FilamentSpecies::FinalizeLocalOrderAnalysis() {
  //if (params_->local_structure_average) {
    //WriteLocalOrderData();
  //}
//}

void FilamentSpecies::FinalizeMse2eAnalysis() {
  int num = members_.size();
  mse2e_file_ << num << " ";
  mse2e_ /= n_samples_*num;
  mse2e2_ /= n_samples_*num;
  mse2e_file_ << mse2e_ << " ";
  mse2e_file_ << sqrt((mse2e2_ - mse2e_*mse2e_)/(num*n_samples_)) << "\n";
}

void FilamentSpecies::FinalizeThetaAnalysis() {
  int nbonds = members_[members_.size()-1].GetNBonds();
  for (int i=0; i<n_bins_; ++i) {
    double axis = (2.0/n_bins_)*i - 1;
    theta_file_ << " " << axis;
    for (int ibond=0; ibond<nbonds-1; ++ibond) {
      theta_file_ << " " << theta_histogram_[ibond][i];
    }
    theta_file_ << "\n";
  }

  for (int ibond=0; ibond<nbonds-1; ++ibond) {
    delete[] theta_histogram_[ibond];
  }
  delete[] theta_histogram_;
}

void FilamentSpecies::InitFlockingAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_filament.flock");
  flock_file_.open(fname, std::ios::out);
  flock_file_ << "flock_analysis_file\n";
  flock_file_ << "length diameter bond_length persistence_length umax driving nsteps nspec delta\n";
  auto it=members_.begin(); 
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  flock_file_ << l << " " << d << " " << cl << " " << pl << " " << params_->soft_potential_mag <<  " " << dr << " " << params_->n_steps << " " << nspec << " " << params_->delta << "\n";
  flock_file_ << "time n_flocking n_interior n_exterior n_joined n_left\n";
}

void FilamentSpecies::RunFlockingAnalysis() {
  int n_flocking = 0;
  int n_interior = 0;
  int n_exterior = 0;
  int n_joined = 0;
  int n_left = 0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    it->CheckFlocking();
    int flock_type = it->GetFlockType();
    int flock_change_state = it->GetFlockChangeState();
    if (flock_type) {
      n_flocking++;
      if (flock_type == 1) {
        n_interior++;
      }
      else if (flock_type == 2) {
        n_exterior++;
      }
      else {
        warning("Flock type not recognized in FilamentSpecies::RunFlockingAnalysis");
      }
    }
    if (flock_change_state) {
      if (flock_change_state == 1) {
        n_joined++;
      }
      else if (flock_change_state == 2) {
        n_left++;
      }
      else {
        warning("Flock change state not recognized in FilamentSpecies::RunFlockingAnalysis");
      }
    }
  }
  if (flock_file_.is_open()) {
    flock_file_ << time_ << " " << n_flocking << " " << n_interior << " " << n_exterior << " " << n_joined << " " << n_left << "\n";
  }
  else {
    early_exit = true;
    std::cout << "ERROR: Problem opening file in RunFlockingAnalysis! Exiting.\n";
  }
}

void FilamentSpecies::FinalizeFlockingAnalysis() {

}
