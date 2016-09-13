// Implementation for simple binding and unbinding

#include <cassert>

#include "xlink_kmc.h"

#include "xlink.h"
#include "xlink_head.h"
#include "xlink_helpers.h"
#include "br_rod.h"

#include <iomanip>

void XlinkKMC::Init(space_struct *pSpace,
                        ParticleTracking *pTracking,
                        PotentialManager *pPotentials,
                        SpeciesBase *spec1,
                        SpeciesBase *spec2,
                        int ikmc,
                        YAML::Node &node,
                        long seed) {
  KMCBase::Init(pSpace, pTracking, pPotentials, spec1, spec2, ikmc, node, seed);

  // Grab our specific claims
  // eps effective
  switch (node["kmc"][ikmc]["concentration_0_1"].Type()) {
    case YAML::NodeType::Scalar:
      eps_eff_0_1_[0] = eps_eff_0_1_[1] =
        0.5 * node["kmc"][ikmc]["concentration_0_1"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      eps_eff_0_1_[0] = node["kmc"][ikmc]["concentration_0_1"][0].as<double>();
      eps_eff_0_1_[1] = node["kmc"][ikmc]["concentration_0_1"][1].as<double>();
      break;
  }
  switch (node["kmc"][ikmc]["concentration_1_2"].Type()) {
    case YAML::NodeType::Scalar:
      eps_eff_1_2_[0] = eps_eff_1_2_[1] =
        0.5 * node["kmc"][ikmc]["concentration_1_2"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      eps_eff_1_2_[0] = node["kmc"][ikmc]["concentration_1_2"][0].as<double>();
      eps_eff_1_2_[1] = node["kmc"][ikmc]["concentration_1_2"][1].as<double>();
      break;
  }

  // on rates
  switch (node["kmc"][ikmc]["on_rate_0_1"].Type()) {
    case YAML::NodeType::Scalar:
      on_rate_0_1_[0] = on_rate_0_1_[1] =
        node["kmc"][ikmc]["on_rate_0_1"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      on_rate_0_1_[0] = node["kmc"][ikmc]["on_rate_0_1"][0].as<double>();
      on_rate_0_1_[1] = node["kmc"][ikmc]["on_rate_0_1"][1].as<double>();
      break;
  }
  switch (node["kmc"][ikmc]["on_rate_1_2"].Type()) {
    case YAML::NodeType::Scalar:
      on_rate_1_2_[0] = on_rate_1_2_[1] =
        node["kmc"][ikmc]["on_rate_1_2"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      on_rate_1_2_[0] = node["kmc"][ikmc]["on_rate_1_2"][0].as<double>();
      on_rate_1_2_[1] = node["kmc"][ikmc]["on_rate_1_2"][1].as<double>();
      break;
  }

  // end pause
  switch (node["kmc"][ikmc]["end_pause"].Type()) {
    case YAML::NodeType::Scalar:
      end_pause_[0] = end_pause_[1] =
        node["kmc"][ikmc]["end_pause"].as<bool>();
      break;
    case YAML::NodeType::Sequence:
      end_pause_[0] = node["kmc"][ikmc]["end_pause"][0].as<bool>();
      end_pause_[1] = node["kmc"][ikmc]["end_pause"][1].as<bool>();
      break;
  }

  // stall force
  switch (node["kmc"][ikmc]["stall_force"].Type()) {
    case YAML::NodeType::Scalar:
      f_stall_[0] = f_stall_[1] =
        node["kmc"][ikmc]["stall_force"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      f_stall_[0] = node["kmc"][ikmc]["stall_force"][0].as<double>();
      f_stall_[1] = node["kmc"][ikmc]["stall_force"][1].as<double>();
      break;
  }

  // velocity
  switch (node["kmc"][ikmc]["velocity"].Type()) {
    case YAML::NodeType::Scalar:
      velocity_[0] = velocity_[1] =
        node["kmc"][ikmc]["velocity"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      velocity_[0] = node["kmc"][ikmc]["velocity"][0].as<double>();
      velocity_[1] = node["kmc"][ikmc]["velocity"][1].as<double>();
      break;
  }

  // velocity polar scale
  switch (node["kmc"][ikmc]["velocity_polar_scale"].Type()) {
    case YAML::NodeType::Scalar:
      velocity_p_scale_[0] = velocity_p_scale_[1] =
        node["kmc"][ikmc]["velocity_polar_scale"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      velocity_p_scale_[0] = node["kmc"][ikmc]["velocity_polar_scale"][0].as<double>();
      velocity_p_scale_[1] = node["kmc"][ikmc]["velocity_polar_scale"][1].as<double>();
      break;
  }

  // velocity antipolar scale
  switch (node["kmc"][ikmc]["velocity_antipolar_scale"].Type()) {
    case YAML::NodeType::Scalar:
      velocity_ap_scale_[0] = velocity_ap_scale_[1] =
        node["kmc"][ikmc]["velocity_antipolar_scale"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      velocity_ap_scale_[0] = node["kmc"][ikmc]["velocity_antipolar_scale"][0].as<double>();
      velocity_ap_scale_[1] = node["kmc"][ikmc]["velocity_antipolar_scale"][1].as<double>();
      break;
  }

  // velocity switch costheta
  switch (node["kmc"][ikmc]["velocity_switch_costheta"].Type()) {
    case YAML::NodeType::Scalar:
      velocity_switch_costheta_[0] = velocity_switch_costheta_[1] =
        node["kmc"][ikmc]["velocity_switch_costheta"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      velocity_switch_costheta_[0] = node["kmc"][ikmc]["velocity_switch_costheta"][0].as<double>();
      velocity_switch_costheta_[1] = node["kmc"][ikmc]["velocity_switch_costheta"][1].as<double>();
      break;
  }

  // Diffusion singly bound
  switch (node["kmc"][ikmc]["diffusion_singly_bound"].Type()) {
    case YAML::NodeType::Scalar:
      diffusion_bound_1_[0] = diffusion_bound_1_[1] =
        node["kmc"][ikmc]["diffusion_singly_bound"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      diffusion_bound_1_[0] = node["kmc"][ikmc]["diffusion_singly_bound"][0].as<double>();
      diffusion_bound_1_[1] = node["kmc"][ikmc]["diffusion_singly_bound"][1].as<double>();
      break;
  }

  // Diffusion doubly bound
  switch (node["kmc"][ikmc]["diffusion_doubly_bound"].Type()) {
    case YAML::NodeType::Scalar:
      diffusion_bound_2_[0] = diffusion_bound_2_[1] =
        node["kmc"][ikmc]["diffusion_doubly_bound"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      diffusion_bound_2_[0] = node["kmc"][ikmc]["diffusion_doubly_bound"][0].as<double>();
      diffusion_bound_2_[1] = node["kmc"][ikmc]["diffusion_doubly_bound"][1].as<double>();
      break;
  }

  alpha_          = node["kmc"][ikmc]["alpha"].as<double>();
  rcutoff_0_1_    = node["kmc"][ikmc]["rcut"].as<double>();
  barrier_weight_ = node["kmc"][ikmc]["barrier_weight"].as<double>();
  k_stretch_      = node["kmc"][ikmc]["spring_constant"].as<double>();
  r_equil_        = node["kmc"][ikmc]["equilibrium_length"].as<double>();
  polar_affinity_ = node["kmc"][ikmc]["polar_affinity"].as<double>();
  write_event_    = node["kmc"][ikmc]["write_event"].as<bool>();

  // XXX Check k_stretch and r_equil
  // done the first time we call the potential

  // Stall type
  std::string stall_str = node["kmc"][ikmc]["stall_type"].as<std::string>();
  if (!stall_str.compare("none")) {
    stall_type_ = 0;
  } else if (!stall_str.compare("parallel")) {
    stall_type_ = 1;
  } else if (!stall_str.compare("absolute")) {
    stall_type_ = 2;
  } else {
    std::cout << "Incorrect stall type: " << stall_str << ", exiting\n";
    exit(1);
  }

  // Things that we need for the tables and CalcCutoff
  BrRodSpecies *prspec = dynamic_cast<BrRodSpecies*>(spec2_);
  max_length_ = prspec->GetMaxLength();

  CalcCutoff();

  BuildTables();
}

void XlinkKMC::CalcCutoff() {
  rcutoff_1_2_ = 0.0;
  const double temp = 1.0;
  const double smalleps = 1E-3;
  double eps_eff = eps_eff_1_2_[0] + eps_eff_1_2_[1];
  double rc_0 = sqrt(2.0 / ( (1-barrier_weight_) * k_stretch_) *
                     temp *
                     log(eps_eff * max_length_ / smalleps * sqrt(2.0 * temp / k_stretch_)));
  rcutoff_1_2_ = r_equil_ + rc_0;
}

double XlinkKMC::XKMCErfinv(double x) {
  // See: A handy approximation for the error function and its inverse
  // (Winitzki 2008) (google it).  This isn't a great approximation, but
  // it will do the trick.  It's not programmed for efficiency since it should
  // only be called during program startup once. If you need better prevision
  // or performance, boost apparently has a version (but then we need boot)
  const double a = 0.147;
  double t1 = -2/M_PI/a;
  double t2 = -log(1-x*x)/2;
  double t3 = 2/M_PI/a + log(1-x*x)/2;
  double t4 = -log(1-x*x)/a;

  return sqrt(t1 + t2 + sqrt(t3*t3 + t4));
}

void XlinkKMC::BuildTables() {
  if (r_equil_ != 0.0) {
    std::vector<double> x[2];
    double bin_size = 0.05;
    double alpha = k_stretch_ * (1 - barrier_weight_) / 2;
    double const smalleps = 1E-5;
    double a_cutoff = 1/sqrt(alpha) * XKMCErfinv(1 - 4.0*sqrt(alpha/M_PI)*smalleps) +
      r_equil_;
    double y_cutoff = rcutoff_1_2_;
    
    xlh::xlink_params params;
    params.alpha = alpha;
    params.r0 = r_equil_;

    for (double a = 0.0; a <= a_cutoff; a += bin_size)
      x[0].push_back(a);
    for (double y0 = 0.0; y0 <= y_cutoff; y0 += bin_size)
      x[1].push_back(y0);

    n_exp_lookup_.Init(2, x, &xlh::prob_1_2, &params);
  }
}

void XlinkKMC::Print() {
  std::cout << "Xlink - BrRod KMC Module\n";
  KMCBase::Print();
  std::cout << std::setprecision(16) << "\teps_eff 0 -> 1:           [" << eps_eff_0_1_[0] << ", " << eps_eff_0_1_[1] << "]\n";
  std::cout << std::setprecision(16) << "\teps_eff 1 -> 2:           [" << eps_eff_1_2_[0] << ", " << eps_eff_1_2_[1] << "]\n";
  std::cout << std::setprecision(16) << "\ton_rate 0 -> 1:           [" << on_rate_0_1_[0] << ", " << on_rate_0_1_[1] << "]\n";
  std::cout << std::setprecision(16) << "\ton_rate 1 -> 2:           [" << on_rate_1_2_[0] << ", " << on_rate_1_2_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tstall_force:              [" << f_stall_[0] << ", " << f_stall_[1] << "]\n";
  std::string stall_str;
  if (stall_type_ == 0) {
    stall_str = "none";
  } else if (stall_type_ == 1) {
    stall_str = "parallel";
  } else if (stall_type_ == 2) {
    stall_str = "absolute";
  } else {
    stall_str = "wtfmate";
  }
  std::cout <<                          "\tstall_type:                " << stall_str << std::endl;
  std::cout <<                          "\tend_pause:                [" << (end_pause_[0] ? "true" : "false") << ", "
                                        << (end_pause_[1] ? "true" : "false") << "]\n";
  std::cout << std::setprecision(16) << "\tvelocity:                 [" << velocity_[0] << ", " << velocity_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tvelocity_polar_scale:     [" << velocity_p_scale_[0] << ", " << velocity_p_scale_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tvelocity_antipolar_scale: [" << velocity_ap_scale_[0] << ", " << velocity_ap_scale_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tvelocity_switch_costheta: [" << velocity_switch_costheta_[0] << ", " << velocity_switch_costheta_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tdiffusion_singly_bound:   [" << diffusion_bound_1_[0] << ", " << diffusion_bound_1_[1] << "]\n";
  std::cout << std::setprecision(16) << "\tdiffusion_doubly_bound:   [" << diffusion_bound_2_[0] << ", " << diffusion_bound_2_[1] << "]\n";
  std::cout << "\tbarrier_weight: " << std::setprecision(16) << barrier_weight_ << std::endl;
  std::cout << "\tequilibrium_length: " << std::setprecision(16) << r_equil_ << std::endl;
  std::cout << "\tk_spring: " << std::setprecision(16) << k_stretch_ << std::endl;
  std::cout << "\tpolar_affinity: " << std::setprecision(16) << polar_affinity_ << std::endl;
  std::cout << "\tmax_length: " << std::setprecision(16) << max_length_ << std::endl;
  std::cout << "\trcutoff_0_1: " << std::setprecision(16) << rcutoff_0_1_ << std::endl;
  std::cout << "\trcutoff_1_2: " << std::setprecision(16) << rcutoff_1_2_ << std::endl;
  std::cout << "\talpha: " << std::setprecision(16) << alpha_ << std::endl;
}

void XlinkKMC::PrepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Prepare each composite particle for the upcoming kmc step
  //if (spec1_->GetSID() != sid1_) return;
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  double ntot_0_1 = 0.0;
  double ntot_1_2 = 0.0;
  auto xlinks = pxspec->GetXlinks(); 

  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Determine what state we're in so that we can do the appropriate thing
    // XXX Do the rest of this
    switch((*xit)->GetBoundState()) {
      case unbound:
        Update_0_1(*xit);
        ntot_0_1 += (*xit)->GetNExp_0_1();
        break;
      case singly:
        Update_1_2(*xit);
        ntot_1_2 += (*xit)->GetNExp_1_2();
        break;
    }
  }

  pxspec->SetNExp_0_1(ntot_0_1);
  pxspec->SetNExp_1_2(ntot_1_2);
}

void XlinkKMC::Update_0_1(Xlink* xit) {
  double nexp_xlink = 0.0;
  auto heads = xit->GetHeads();

  // Final binding affinity is 2x if the eps_eff and on_rate are the same
  for (int i = 0; i < heads->size(); ++i) {
    auto head = &(*heads)[i];
    double nexp = 0.0;
    double binding_affinity = eps_eff_0_1_[i] * on_rate_0_1_[i] * alpha_ * xit->GetDelta();
    auto idx = (*oid_position_map_)[head->GetOID()];
    for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
      nexp += binding_affinity * nldx->kmc_;
    }
    head->SetNExp_0_1(nexp);
    nexp_xlink += nexp;
  }

  xit->SetNExp_0_1(nexp_xlink);
}

void XlinkKMC::Update_1_2(Xlink *xit) {
  XlinkHead *freehead, *boundhead;
  auto isbound = xit->GetBoundHeads(&freehead, &boundhead);
  // If the first head is bound, then the second one is free, and vice
  // versa
  int free_i = isbound.first ? 1 : 0;

  double binding_affinity = eps_eff_1_2_[free_i] * on_rate_1_2_[free_i];
  auto free_idx = (*oid_position_map_)[freehead->GetOID()];
  // Get the attached rod
  auto attach_info = boundhead->GetAttach();
  auto attach_info_idx = (*oid_position_map_)[attach_info.first];
  auto mrod_attached = (*simples_)[attach_info_idx];
  if (binding_affinity > 0.0) {
    double n_exp = 0.0;

    // We have to look at all of our neighbors withint the mrcut
    for (auto nldx = neighbors_[free_idx].begin(); nldx != neighbors_[free_idx].end(); ++nldx) {
      auto mrod = (*simples_)[nldx->idx_];
      // Check to see if it's really a rod, and if it's the same one we're already attached to
      if (mrod->GetSID() != sid2_) continue;
      if (mrod->GetRID() == mrod_attached->GetRID()) {
        continue;
      }

      // Calculate center to center displacement
      double polar_affinity = xlh::polar_affinity(ndim_, polar_affinity_, mrod_attached->GetRigidOrientation(), mrod->GetRigidOrientation());
      double r_x[3];
      double s_x[3];
      double r_rod[3];
      double s_rod[3];
      double u_rod[3];
      std::copy(freehead->GetRigidPosition(), freehead->GetRigidPosition()+ndim_, r_x);
      std::copy(freehead->GetRigidScaledPosition(), freehead->GetRigidScaledPosition()+ndim_, s_x);
      std::copy(mrod->GetRigidPosition(), mrod->GetRigidPosition()+ndim_, r_rod);
      std::copy(mrod->GetRigidScaledPosition(), mrod->GetRigidScaledPosition()+ndim_, s_rod);
      std::copy(mrod->GetRigidOrientation(), mrod->GetRigidOrientation()+ndim_, u_rod);
      double l_rod = mrod->GetRigidLength();
      double dr[3] = {0.0, 0.0, 0.0};
      double mu0 = 0.0;

      min_distance_point_carrier_line_inf(ndim_, nperiodic_,
                                          space_->unit_cell,
                                          r_x, s_x,
                                          r_rod, s_rod, u_rod, l_rod,
                                          dr, &mu0);

      // Now do the integration over the limits on the MT
      double r_min_mag2 = 0.0;
      // Calculated across the space between the min distance and the opposing
      for (int i = 0; i < ndim_; ++i) {
        double dri = u_rod[i] * mu0 + dr[i];
        r_min_mag2 += SQR(dri);
      }
      // Check cutoff distance
      if (r_min_mag2 > SQR(rcutoff_1_2_)) {
        nldx->kmc_ = 0.0;
        n_exp += nldx->kmc_;
        continue;
      }
      if (r_equil_ == 0.0) {
        double kb = k_stretch_ * (1.0 - barrier_weight_);
        double scale_factor = sqrt(0.5 * kb);
        double lim0 = scale_factor * (-mu0 - 0.5*l_rod);
        double term0 = erf(lim0);
        double lim1 = scale_factor * (-mu0 + 0.5*l_rod);
        double term1 = erf(lim1);

        nldx->kmc_ = binding_affinity * sqrt(M_PI_2 / kb) * exp(-0.5*kb*r_min_mag2) * (term1 - term0) * polar_affinity;
        n_exp += nldx->kmc_;
      } else {
        double lim0 = -mu0 - 0.5 * l_rod;
        double lim1 = -mu0 + 0.5 * l_rod;
        double r_min_mag = sqrt(r_min_mag2);
        double x[2] = {fabs(lim0), r_min_mag};
        double term0 = n_exp_lookup_.Lookup(x) * ((lim0 < 0) ? -1.0 : 1.0);
        x[0] = fabs(lim1);
        double term1 = n_exp_lookup_.Lookup(x) * ((lim1 < 0) ? -1.0 : 1.0);
        // OVERRIDE the kmc_ value of this neighbor list
        nldx->kmc_ = binding_affinity * (term1 - term0) * polar_affinity;
        n_exp += nldx->kmc_;
      }
      if (debug_trace) {
        std::cout << "[" << freehead->GetOID() << "] -> neighbor[" << mrod->GetOID() << "]";
        std::cout << " {kmc: " << std::setprecision(16) << nldx->kmc_ << "}\n";
      }
    } // loop over local neighbors of xlink

    boundhead->SetNExp_1_2(0.0);
    freehead->SetNExp_1_2(n_exp);
    xit->SetNExp_1_2(n_exp);
  }
}

void XlinkKMC::StepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Run the bind unbind
  int g[4] = {0, 1, 2, 3};
  for (int i = 0; i < 4; ++i) {
    int j = (int)gsl_rng_uniform_int(rng_.r, 4);
    int swapme = g[i];
    g[i] = g[j];
    g[j] = swapme;
  }

  if (debug_trace)
    std::cout << "XlinkKMC module " << g[0] << " -> " << g[1] << " -> " << g[2] << " -> " << g[3] << std::endl;

  for (int i = 0; i < 4; ++i) {
    switch (g[i]) {
      case 0:
        KMC_0_1();
        break;
      case 1:
        KMC_1_0();
        break;
      case 2:
        KMC_1_2();
        break;
      case 3:
        KMC_2_1();
        break;
    }
  }
}

void XlinkKMC::KMC_0_1() {
  // Loop over the xlinks to see who binds
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Only take free ones
    if ((*xit)->GetBoundState() != unbound) continue;
    auto nexp = (*xit)->GetNExp_0_1(); // Number expected to bind in this timestep, calculated from Update_0_1
    if (nexp <  std::numeric_limits<double>::epsilon() &&
        nexp > -std::numeric_limits<double>::epsilon()) nexp = 0.0;
    // IF we have some probability to fall onto a neighbor, check it
    if (nexp > 0.0) {
      std::cout << "--------\n";
      std::cout << "Xlink: " << (*xit)->GetOID() << std::endl;
      auto mrng = (*xit)->GetRNG();
      double roll = gsl_rng_uniform(mrng->r);
      std::cout << std::setprecision(16) << "roll: " << roll << std::endl;
      if (roll < nexp) {
        std::cout << "Successful roll\n";
        int head_type = gsl_rng_uniform(mrng->r) < ((eps_eff_0_1_[1])/(eps_eff_0_1_[0]+eps_eff_0_1_[1]));
        std::cout << "head: " << head_type << std::endl;
        auto heads = (*xit)->GetHeads();
        auto head = heads->begin() + head_type;
        double binding_affinity = (eps_eff_0_1_[0] * on_rate_0_1_[0] + eps_eff_0_1_[1] * on_rate_0_1_[1]) *
          alpha_ * head->GetDelta();
        std::cout << "binding_affinity: " << std::setprecision(16) << binding_affinity << std::endl;
        std::ostringstream kmc_event;
        kmc_event << std::setprecision(16) << "[" << (*xit)->GetOID() << "] Successful KMC move {0 -> 1}, {nexp: " << nexp;
        kmc_event << std::setprecision(16) << "}, {roll: " << roll << "}, {head: " << head_type << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace) {
          std::cout << kmc_event.str() << std::endl;
        }
        double pos = 0.0;
        int idx = (*oid_position_map_)[head->GetOID()];
        // Search through the neighbors of this head to figure out who we want to bind to
        int ineighb = 0;
        // DEBUG XXX FIXME totalval
        double totalval = 0.0;
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
          totalval += binding_affinity * nldx->kmc_;
        }
        std::cout << "totalval: " << std::setprecision(16) << totalval << std::endl;
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx, ++ineighb) {
          std::cout << "[" << idx << "] -> neighbor[" << ineighb << "]\n";
          auto part2 = (*simples_)[nldx->idx_];
          if (part2->GetSID() != sid2_) continue; // Make sure it's what we want to bind to
          std::cout << std::setprecision(16) << "my binding aff: " << binding_affinity * nldx->kmc_ << std::endl;
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            if (debug_trace) {
              std::cout << std::setprecision(16) << "[" << head->GetOID() << "] Attaching to ["
                << part2->GetOID() << "], Z: " << std::setprecision(16) << pos << std::endl;
            }

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double s_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(head->GetRigidPosition(), head->GetRigidPosition()+ndim_, r_x);
            std::copy(head->GetRigidScaledPosition(), head->GetRigidScaledPosition()+ndim_, s_x);
            std::copy(part2->GetRigidPosition(), part2->GetRigidPosition()+ndim_, r_rod);
            std::copy(part2->GetRigidScaledPosition(), part2->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(part2->GetRigidOrientation(), part2->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = part2->GetRigidLength();
            double dr[3];
            double mu = 0.0;
            min_distance_point_carrier_line_inf(ndim_, nperiodic_,
                                                space_->unit_cell, r_x, s_x,
                                                r_rod, s_rod, u_rod, l_rod,
                                                dr, &mu);

            double r_min[3];
            double r_min_mag2 = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              r_min[i] = -mu * u_rod[i] - dr[i];
              r_min_mag2 += SQR(r_min[i]);
            }
            mrcut2_ = rcutoff_0_1_*rcutoff_0_1_;
            double a = sqrt(mrcut2_ - r_min_mag2);
            //double a = sqrt(1.0 - r_min_mag2); //FIXME is this right for 1.0? or mrcut2?
            if (a!=a)
              a = 0.0;

            std::cout << std::setprecision(16) << "rminmag2: " << r_min_mag2 << ", a: " << a << std::endl;

            double crosspos = 0.0;
            for (int i = 0; i < 100; ++i) {
              double uroll = gsl_rng_uniform(mrng->r);
              crosspos = (uroll - 0.5)*a + mu + 0.5*l_rod;

              if (crosspos >= 0 && crosspos <= l_rod)
                break;
              if (i == 99) {
                crosspos = -mu + 0.5*l_rod;
                if (crosspos < 0)
                  crosspos = 0;
                else if (crosspos > l_rod)
                  crosspos = l_rod;
              }
            }
            // Attach to the OID of the particle, this is done for dynamic instability to work
            head->Attach(part2->GetOID(), crosspos);
            if (debug_trace) {
              std::cout << std::setprecision(16) << "\t{" << mu << "}, {crosspos: " << crosspos << "}\n";
            }
            head->SetBound(true);
            (*xit)->CheckBoundState();
            break;
          } // the actual neighbor to fall on
        } // loop over neighbors to figure out which to attach to
      } // successfully bound
    } // found an xlink with some expectation of binding
  } // loop over all xlinks
}

void XlinkKMC::KMC_1_0() {
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  int nbound1[2];
  std::copy(pxspec->GetNBound1(), pxspec->GetNBound1()+2, nbound1);
  double poff_single_a = on_rate_0_1_[0] * alpha_ * pxspec->GetDelta();
  double poff_single_b = on_rate_0_1_[1] * alpha_ * pxspec->GetDelta();

  std::cout << std::setprecision(16) << "poff: [" << poff_single_a << ", " << poff_single_b << "]\n";
  std::cout << "nbound1: [" << nbound1[0] << ", " << nbound1[1] << "]\n";

  int noff[2] = {(int)gsl_ran_binomial(rng_.r, poff_single_a, nbound1[0]),
                 (int)gsl_ran_binomial(rng_.r, poff_single_b, nbound1[1])};
  if (debug_trace) {
    std::cout << std::setprecision(16) << "[KMC_1_0] poff_single: [" << poff_single_a << ", " << poff_single_b
      << "], Removing [" << noff[0] << ", " << noff[1] << "]\n";
  }

  for (int i = 0; i < (noff[0] + noff[1]); ++i) {
    if (debug_trace) {
      std::cout << "[KMC_1_0] detaching trial " << i+1 << "/" << (noff[0] + noff[1]) << std::endl;
    }
    int head_type = i < noff[1];
    int idxloc = -1;
    int idxoff = (int)gsl_rng_uniform_int(rng_.r, nbound1[head_type]);

    std::cout << "head_type: " << head_type << ", idxoff " << idxoff << std::endl;

    // Find the one to remove
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      if ((*xit)->GetBoundState() != singly) continue;
      auto heads = (*xit)->GetHeads();
      // Check to increment only in the case we have the correctly
      // bound head
      if (!(*heads)[head_type].GetBound()) continue;
      idxloc++;
      if (idxloc == idxoff) {
        // If we've found it, make sure that we decrement nbound1
        // so that we get the correct counting in the case that we remove more
        // than 1 head in a single timestep
        //
        // For example, we could remove the last idx, and then remove the last
        // idx again (With some change), and we want to make sure we actually
        // do detach the correct number
        nbound1[head_type]--;
        XlinkHead *boundhead, *freehead;
        auto isbound = (*xit)->GetBoundHeads(&freehead, &boundhead);

        std::ostringstream kmc_event;
        kmc_event << "[x:" << (*xit)->GetOID() << ",head:" << boundhead->GetOID() << "] Successful KMC move {1 -> 0}, {idxoff=idxloc=";
        kmc_event << idxloc << "}, {head: " << head_type << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace) {
          std::cout << kmc_event.str() << std::endl;
        }

        // Call the single head detach function
        Detach_1_0((*xit), freehead, boundhead);
        break;

      } // found it!
    } // find the one to remove

  } // How many to remove
}

void XlinkKMC::KMC_1_2() {
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  double nexp_1_2 = pxspec->GetNExp_1_2();
  if (debug_trace)
    printf("[KMC_1_2] ntot: %2.8f\n", nexp_1_2 * pxspec->GetDelta());
  int nattach = gsl_ran_poisson(rng_.r, nexp_1_2 * pxspec->GetDelta());
  for (int itrial = 0; itrial < nattach; ++itrial) {
    if (debug_trace)
      printf("[KMC_1_2] attaching trial %d/%d\n", itrial+1, nattach);
    double ran_loc = gsl_rng_uniform(rng_.r) * nexp_1_2;
    double loc = 0.0;
    bool foundidx = false;
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      // Only take singly bound
      if ((*xit)->GetBoundState() != singly) continue;
      loc += (*xit)->GetNExp_1_2();
      if (loc > ran_loc) {
        std::ostringstream kmc_event;
        kmc_event << "[" << (*xit)->GetOID() << "] Successful KMC move {1 -> 2}, {nexp_1_2: " << nexp_1_2 << "}, {ran_loc: ";
        kmc_event << ran_loc << "}, {loc: " << loc << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace)
          printf("%s\n", kmc_event.str().c_str());
        // reset loc back to beginning of this xlink
        loc -= (*xit)->GetNExp_1_2();

        // Look at my neighbors and figure out which to fall on
        // Also, get the heads
        XlinkHead *freehead, *boundhead;
        auto isbound = (*xit)->GetBoundHeads(&freehead, &boundhead);
        int idx = (*oid_position_map_)[freehead->GetOID()];

        // Loop over my neighbors to figure out which to fall on
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
          auto mrod = (*simples_)[nldx->idx_];
          // Make sure it's the right species to attach to
          if (mrod->GetSID() != sid2_) continue;

          // Ignore our own rod
          //auto free_idx = (*oid_position_map_)[nonattachead->GetOID()];
          //auto free_idx = idx;
          //auto attc_idx = (*oid_position_map_)[attachedhead->GetOID()];
          // Get the attached rod
          auto attachinfo = boundhead->GetAttach();
          auto attachinfoidx = (*oid_position_map_)[attachinfo.first];
          auto mrod_attached = (*simples_)[attachinfoidx];
          if (mrod->GetRID() == mrod_attached->GetRID()) {
            continue;
          }

          loc += nldx->kmc_;
          if (loc > ran_loc) {
            // Found the neighbor to fall onto!!!!
            if (debug_trace) {
              printf("[%d] Attaching to [%d] {loc: %2.4f}\n", freehead->GetOID(), mrod->GetOID(), loc);
            }

            // What location (crosspos) do we fall onto?
            double crosspos = 0.0;

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double s_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(freehead->GetRigidPosition(), freehead->GetRigidPosition()+ndim_, r_x);
            std::copy(freehead->GetRigidScaledPosition(), freehead->GetRigidScaledPosition()+ndim_, s_x);
            std::copy(mrod->GetRigidPosition(), mrod->GetRigidPosition()+ndim_, r_rod);
            std::copy(mrod->GetRigidScaledPosition(), mrod->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(mrod->GetRigidOrientation(), mrod->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = mrod->GetRigidLength();
            double rcontact[3];
            double dr[3];
            double mu = 0.0;
            min_distance_point_carrier_line(ndim_, nperiodic_,
                                            space_->unit_cell, r_x, s_x,
                                            r_rod, s_rod, u_rod, l_rod,
                                            dr, rcontact, &mu);

            // Find it
            auto mrng = (*xit)->GetRNG();
            if (r_equil_ == 0.0) {
              double kb = (1.0 - barrier_weight_) * k_stretch_;
              do {
                double mpos = gsl_ran_gaussian_ziggurat(mrng->r, sqrt(1.0/kb)) + mu + 0.5 * l_rod;
                if (mpos >= 0 && mpos <= l_rod) {
                  crosspos = mpos;
                  break;
                }
              } while(1);
            } else {
              double y02 = 0.0;
              for (int i = 0; i < ndim_; ++i) {
                y02 += SQR(dr[i]);
              }
              int itrial_loc = 0;
              do {
                itrial_loc++;
                double uroll = gsl_rng_uniform(mrng->r);
                double xvec[2] = {0.0, sqrt(y02)};
                double mpos = ((gsl_rng_uniform(mrng->r) < 0.5) ? -1.0 : 1.0) *
                    n_exp_lookup_.Invert(0, uroll, xvec) + mu + 0.5 * l_rod;
                if (mpos >= 0 && mpos <= l_rod) {
                  crosspos = mpos;
                  break;
                }
              } while (itrial_loc < 100);
            }

            freehead->Attach(mrod->GetOID(), crosspos);
            freehead->SetBound(true);
            (*xit)->CheckBoundState();

            // Update the position of the xlink for the force calculation
            auto headid_free = freehead->GetHeadID();
            auto headid_bound = boundhead->GetHeadID();

            if (headid_free == 0 && headid_bound == 1) {
              (*xit)->UpdateStagePosition(mrod->GetRigidPosition(),
                                          mrod->GetRigidOrientation(),
                                          mrod->GetRigidLength(),
                                          mrod->GetOID(),
                                          mrod_attached->GetRigidPosition(),
                                          mrod_attached->GetRigidOrientation(),
                                          mrod_attached->GetRigidLength(),
                                          mrod_attached->GetOID());
            } else if (headid_free == 1 && headid_bound == 0) {
              (*xit)->UpdateStagePosition(mrod_attached->GetRigidPosition(),
                                          mrod_attached->GetRigidOrientation(),
                                          mrod_attached->GetRigidLength(),
                                          mrod_attached->GetOID(),
                                          mrod->GetRigidPosition(),
                                          mrod->GetRigidOrientation(),
                                          mrod->GetRigidLength(),
                                          mrod->GetOID());
            } else {
              std::cout << "Some attachment has gone horribly wrong!\n";
              exit(1);
            }

            // Calculate the potentials and forces of this xlink
            PotentialBase *xlink_pot = potentials_->GetPotentialInternal(freehead->GetOID(), boundhead->GetOID());
            if (xlink_pot == nullptr) {
              std::cout << "Uhhhh......\n";
              exit(1);
            }

            if (!first_potential_use) {
              first_potential_use = true;
              XlinkHarmonic *xharm = dynamic_cast<XlinkHarmonic*>(xlink_pot);
              auto myk = xharm->GetK();
              auto myrequil = xharm->GetRequil();
              if (myk != k_stretch_) {
                std::cout << "Uh oh, k: " << k_stretch_ << " != potential k: " << myk << std::endl;
                exit(1);
              }
              if (myrequil != r_equil_) {
                std::cout << "Uh oh, requil: " << r_equil_ << " != potential requil: " << myrequil << std::endl;
                exit(1);
              }
            }

            // Calculate the minimum distance, regardless of any cutoff
            interactionmindist idm;
            MinimumDistance(freehead, boundhead, idm, ndim_, nperiodic_, space_);

            // Fire off the potential calculation
            double fepot[4];
            xlink_pot->CalcPotential(&idm, freehead, boundhead, fepot);
           
            double flink_free[3] = {0.0, 0.0, 0.0};
            double flink_bound[3] = {0.0, 0.0, 0.0};
            for (int idim = 0; idim < ndim_; ++idim) {
              flink_free[idim] = fepot[idim];
              flink_bound[idim] = -fepot[idim];
            }

            freehead->AddPotential(fepot[ndim_]);
            freehead->AddForce(flink_free);
            boundhead->AddPotential(fepot[ndim_]);
            boundhead->AddForce(flink_bound);

            foundidx = true;
            break;
          } // got the neighbor to fall onto
        } // check neighbors to see if we need to fall onto this one

      } // found it!

      if (foundidx)
        break;
    } // loop over all xlinks for itrial
  }
}

void XlinkKMC::KMC_2_1() {
  if (barrier_weight_ == 0.0) {
    // All detachments equally probable, do via a poisson distribution
    KMC_2_1_ForceIndep();
  } else {
    KMC_2_1_ForceDep();
  }
}

void XlinkKMC::KMC_2_1_ForceIndep() {
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  int nheads[2] = {0, 0};
  nheads[0] = pxspec->GetNBound2()[0];
  nheads[1] = pxspec->GetNBound2()[1];

  int noff[2] = {
    (int)gsl_ran_binomial(rng_.r, on_rate_1_2_[0] * pxspec->GetDelta(), nheads[0]),
    (int)gsl_ran_binomial(rng_.r, on_rate_1_2_[1] * pxspec->GetDelta(), nheads[1])};
  if (debug_trace)
    printf("[Xlink] {poff_head: (%2.8f, %2.8f)}, {noff: (%d, %d)}\n",
        on_rate_1_2_[0] * pxspec->GetDelta(), on_rate_1_2_[1] * pxspec->GetDelta(),
        noff[0], noff[1]);

  for (int itrial = 0; itrial < (noff[0] + noff[1]); ++itrial) {
    if (debug_trace)
      printf("[KMC_2_1] detaching trial %d/%d\n", itrial+1, (noff[0] + noff[1]));
    int head_type = itrial < noff[1];
    int idxloc = -1;
    int idxoff = (int)gsl_rng_uniform_int(rng_.r, nheads[head_type]);
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      if ((*xit)->GetBoundState() != doubly) continue;
      idxloc++;
      if (idxloc == idxoff) {
        nheads[head_type]--;

        std::ostringstream kmc_event;
        kmc_event << "[" << (*xit)->GetOID() << ",head: " << head_type << "] Successful KMC move {2 -> 1}, {head_type: ";
        kmc_event << head_type << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace)
          printf("%s\n", kmc_event.str().c_str());

        // Have to see if both heads detaching, or just one (ugh)
        XlinkHead *freehead, *boundhead;
        auto isbound = (*xit)->GetBoundHeads(&freehead, &boundhead);
        if (isbound.first && isbound.second) {
          Detach_2_1(*xit, head_type);
        } else if (isbound.first && !isbound.second) {
          Detach_1_0(*xit, freehead, boundhead);
        } else if (isbound.second && !isbound.first) {
          Detach_1_0(*xit, freehead, boundhead);
        } else {
          printf("Something has gone horribly wrong!\n");
          exit(1);
        }

        break;
      }
    } // loop over xlinks

  } // remove ntrials
}

void XlinkKMC::KMC_2_1_ForceDep() {
  // Force dependent detachment!
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  double poff_base[2] = {
    on_rate_1_2_[0] * pxspec->GetDelta(),
    on_rate_1_2_[1] * pxspec->GetDelta()};
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    if ((*xit)->GetBoundState() != doubly) continue; // only doubly bound
    double kboltzoff = exp(barrier_weight_ * (*xit)->GetInternalEnergy());
    //printf("kboltzoff: %2.8f\n", kboltzoff);
    // XXX Possibly do this via the xlink rng, multithreaded
    auto xrng = (*xit)->GetRNG();
    uint8_t off[2] = {
      gsl_rng_uniform(xrng->r) < poff_base[0] * kboltzoff,
      gsl_rng_uniform(xrng->r) < poff_base[1] * kboltzoff};
    if (debug_trace)
      printf("[Xlink] {poff_base: (%2.8f, %2.8f)}, {weighted: (%2.8f, %2.8f)}\n",
          poff_base[0], poff_base[1], poff_base[0]*kboltzoff, poff_base[1]*kboltzoff);
    if (off[0] && off[1]) {
      // Call the detach function
      std::ostringstream kmc_event;
      kmc_event << "[" << (*xit)->GetOID() << "] Successful KMC move {2 -> 0}";
      WriteEvent(kmc_event.str());
      if (debug_trace)
        printf("%s\n", kmc_event.str().c_str());

      Detach_2_0(*xit);
    } else if (off[0] || off[1]) {
      // Figure out which head, then call the single detach
      int whichhead = 0;
      if (off[0] == 1) {
        whichhead = 0;
      } else if (off[1] == 1) {
        whichhead = 1;
      } else {
        printf("Something has gone horribly wrong!\n");
        exit(1);
      }

      std::ostringstream kmc_event;
      kmc_event << "[" << (*xit)->GetOID() << ",head: " << whichhead << "] Successful KMC move {2 -> 1}, {off: (";
      kmc_event << (int)off[0] << ", " << (int)off[1] << ")}";
      WriteEvent(kmc_event.str());
      if (debug_trace)
        printf("%s\n", kmc_event.str().c_str());

      Detach_2_1(*xit, whichhead);
    }
  }
}

void XlinkKMC::Detach_1_0(Xlink *xit, XlinkHead *freehead, XlinkHead *boundhead) {
  // Let's get the information and make sure it's right
  std::ostringstream kmc_event;
  kmc_event << "    [" << xit->GetOID() << "] Singly Bound Single Head Detach {" << boundhead->GetOID() << "}";
  WriteEvent(kmc_event.str());
  if (debug_trace) {
    std::cout << kmc_event.str() << std::endl;
  }
  // Place within some random distance of the attach point
  double randr[3];
  double mag2 = 0.0;
  mrcut2_ = rcutoff_0_1_*rcutoff_0_1_;
  auto mrng = boundhead->GetRNG();
  double prevpos[3] = {0.0, 0.0, 0.0};
  std::copy(boundhead->GetRigidPosition(), boundhead->GetRigidPosition()+ndim_, prevpos);
  do {
    mag2 = 0.0;
    for (int i = 0; i < ndim_; ++i) {
      double mrand = gsl_rng_uniform(mrng->r);
      randr[i] = 2*rcutoff_0_1_*(mrand - 0.5);
      mag2 += SQR(randr[i]);
    }
  } while(mag2 > mrcut2_);
  // Randomly set position based on randr
  for (int i = 0; i < ndim_; ++i) {
    randr[i] = randr[i] + prevpos[i];
  }

  boundhead->SetPosition(randr);
  boundhead->SetPrevPosition(prevpos);
  boundhead->UpdatePeriodic();
  boundhead->AddDr();
  freehead->SetPosition(randr);
  freehead->SetPrevPosition(prevpos);
  freehead->UpdatePeriodic();
  freehead->AddDr();

  if (debug_trace) {
    auto attachid = boundhead->GetAttach();
    auto attachidx = (*oid_position_map_)[attachid.first];
    auto part2 = (*simples_)[attachidx];
    std::cout << "[" << xit->GetOID() << "]{" << boundhead->GetOID() << "}";
    std::cout << " Detached from [" << part2->GetOID() << "]";
    std::cout << " (";
    std::cout << std::setprecision(16) << prevpos[0] << ", ";
    std::cout << std::setprecision(16) << prevpos[1] << ", ";
    std::cout << std::setprecision(16) << prevpos[2] << ") ->";
    std::cout << " (";
    std::cout << std::setprecision(16) << boundhead->GetRigidPosition()[0] << ", ";
    std::cout << std::setprecision(16) << boundhead->GetRigidPosition()[1] << ", ";
    std::cout << std::setprecision(16) << boundhead->GetRigidPosition()[2] << ")\n";
  }

  boundhead->SetBound(false);
  boundhead->Attach(-1, 0.0);
  xit->CheckBoundState();
}

void XlinkKMC::Detach_2_1(Xlink *xit, int headtype) {
  // Remove head headtype, set location to other head (not double detach)
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;

  // Figure out which head detaches
  // Do some fancy aliasing to make this easier
  XlinkHead *attachedhead;
  XlinkHead *detachedhead;
  if (headtype == 0) {
    attachedhead = &(*head1);
    detachedhead = &(*head0);
  } else if (headtype == 1) {
    attachedhead = &(*head0);
    detachedhead = &(*head1);
  } else {
    printf("Something has gone horribly wrong!\n");
    exit(1);
  }

  std::ostringstream kmc_event;
  kmc_event << "    [" << xit->GetOID() << "] Doubly Bound Single Head Detach [" << detachedhead->GetOID() << "]";
  kmc_event << " ([" << attachedhead->GetOID() << "] still bound)";
  WriteEvent(kmc_event.str());
  if (debug_trace)
    std::cout << kmc_event.str() << std::endl;

  // Just set the location to what the other head was, and unset the attachment
  double oldpos[3] = {0.0, 0.0, 0.0};
  std::copy(detachedhead->GetRigidPosition(), detachedhead->GetRigidPosition()+ndim_, oldpos);

  detachedhead->SetPosition(attachedhead->GetRigidPosition());
  detachedhead->SetPrevPosition(oldpos);
  detachedhead->UpdatePeriodic();
  detachedhead->AddDr();

  if (debug_trace) {
    auto attachid = detachedhead->GetAttach();
    auto ridx = (*oid_position_map_)[attachid.first];
    auto mrod = (*simples_)[ridx];
   
    std::cout << "[" << xit->GetOID() << "]{" << detachedhead->GetOID() << "}";
    std::cout << " Detached from [" << mrod->GetOID() << "] -> (";
    std::cout << std::setprecision(16) << oldpos[0] << ", " << oldpos[1] << ", " << oldpos[2] << ")"
      << " -> (" << detachedhead->GetRigidPosition()[0] << ", "
      << detachedhead->GetRigidPosition()[1] << ", "
      << detachedhead->GetRigidPosition()[2] << ")\n";
  }

  detachedhead->SetBound(false);
  detachedhead->Attach(-1, 0.0);
  // Update main xlink location
  xit->SetPosition(attachedhead->GetRigidPosition());
  xit->SetPrevPosition(oldpos);
  xit->UpdatePeriodic();
  xit->CheckBoundState();
}

void XlinkKMC::Detach_2_0(Xlink *xit) {
  // Both heads are detaching, set the location to the midpoint of the two heads
  // Easy, set the head to the midpoint of the xlink and fall off both
  auto heads = xit->GetHeads();
  XlinkHead *head0 = &(*(heads->begin()));
  XlinkHead *head1 = &(*(heads->begin()+1));

  std::ostringstream kmc_event;
  kmc_event << "    [" << xit->GetOID() << "] Doubly Bound Double Head Detach [" << head0->GetOID() << ", ";
  kmc_event << head1->GetOID() << "]";
  WriteEvent(kmc_event.str());
  if (debug_trace)
    std::cout << kmc_event.str() << std::endl;

  double oldpos0[3] = {0.0, 0.0, 0.0};
  double oldpos1[3] = {0.0, 0.0, 0.0};
  double xpos[3] = {0.0, 0.0, 0.0};
  std::copy(xit->GetPosition(), xit->GetPosition()+ndim_, xpos);
  std::copy(head0->GetRigidPosition(), head0->GetRigidPosition()+ndim_, oldpos0);
  std::copy(head1->GetRigidPosition(), head1->GetRigidPosition()+ndim_, oldpos1);

  if (debug_trace) {
    auto head0attach = head0->GetAttach();
    auto mrod0idx = (*oid_position_map_)[head0attach.first];
    auto mrod0 = (*simples_)[mrod0idx];
    auto head1attach = head1->GetAttach();
    auto mrod1idx = (*oid_position_map_)[head1attach.first];
    auto mrod1 = (*simples_)[mrod1idx];

    std::cout << "[" << xit->GetOID() << "]{" << head0->GetOID() << "," << head1->GetOID() << "}";
    std::cout << " Detached from [" << mrod0->GetOID() << "," << mrod1->GetOID() << "] -> {";
    std::cout << std::setprecision(16) << "(" << oldpos0[0] << ", " << oldpos0[1] << ", " << oldpos0[2] << "),("
      << oldpos1[0] << ", " << oldpos1[1] << ", " << oldpos1[2] << ")} -> (" << xpos[0] << ", " << xpos[1]
      << ", " << xpos[2] << ")\n";
  }

  head0->SetPosition(xpos);
  head0->SetPrevPosition(oldpos0);
  head0->UpdatePeriodic();
  head0->SetBound(false);
  head0->Attach(-1, 0.0);
  head0->AddDr();
  head1->SetPosition(xpos);
  head1->SetPrevPosition(oldpos1);
  head1->UpdatePeriodic();
  head1->SetBound(false);
  head1->Attach(-1, 0.0);
  head1->AddDr();
  xit->CheckBoundState();
}

void XlinkKMC::UpdateKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  nfree_ = 0.0;
  nbound1_[0] = nbound1_[1] = 0.0;
  nbound2_[0] = nbound2_[1] = 0.0;

  // Do a switch on the type that we're examining
  // Loop over the xlinks to see who  does what
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    switch((*xit)->GetBoundState()) {
      case unbound:
        nfree_++;
        break;
      case singly:
        UpdateStage1(*xit);
        break;
      case doubly:
        UpdateStage2(*xit);
        break;
    }
  }

  pxspec->SetNFree(nfree_);
  pxspec->SetNBound1(nbound1_[0], nbound1_[1]);
  pxspec->SetNBound2(nbound2_[0], nbound2_[0]);
}

void XlinkKMC::UpdateStage1(Xlink *xit) {

  // Set nexp to zero for all involved
  xit->SetNExp_0_1(0.0);
  XlinkHead *freehead, *boundhead;
  auto isbound = xit->GetBoundHeads(&freehead, &boundhead);
  freehead->SetNExp_0_1(0.0);
  boundhead->SetNExp_0_1(0.0);

  bool unbind = false;
  int bound_idx = 0;

  if (isbound.first) {
    bound_idx = 0;
  } else if (isbound.second) {
    bound_idx = 1;
  } else {
    printf("Something has gone horribly wrong!\n");
    exit(1);
  }

  auto aid = boundhead->GetAttach().first;
  auto aidx = (*oid_position_map_)[aid];
  auto cross_pos = boundhead->GetAttach().second; // relative to the -end of the rod!!
  
  auto part2 = (*simples_)[aidx];
  auto l_rod = part2->GetRigidLength();

  // Check the diffusion
  if (diffusion_bound_1_[bound_idx] > 0.0) {
    cross_pos += sqrt(2.0 * diffusion_bound_1_[bound_idx] * boundhead->GetDelta()) *
      gsl_ran_gaussian_ziggurat(boundhead->GetRNG()->r, 1.0);
  }

  // If we are moving with some velocity, do that
  // also check for end pausing
  cross_pos += velocity_[bound_idx] * boundhead->GetDelta();
  if (cross_pos > l_rod) {
    cross_pos = l_rod;
    if (!end_pause_[bound_idx]) {
      unbind = true;
    }
  } else if (cross_pos < 0.0) {
    cross_pos = 0.0;
    if (!end_pause_[bound_idx]) {
      unbind = true;
    }
  }

  // Still need the information on where we were (end of the rod to detach, if possible)
  boundhead->Attach(aid, cross_pos);
  xit->UpdateStagePosition(part2->GetRigidPosition(),
                           part2->GetRigidOrientation(),
                           part2->GetRigidLength(),
                           part2->GetOID(),
                           nullptr,
                           nullptr,
                           0.0,
                           0);
  if (unbind) {
    // Detach this head and put the xlink back out into the nucleoplasm
    Detach_1_0(xit, freehead, boundhead);
  } else {
    // Head is still attached, so we're all good
    nbound1_[bound_idx]++;
  }
}

void XlinkKMC::UpdateStage2(Xlink *xit) {
  // Set nexp to zero for all involved
  // New hotness
  xit->SetNExp_1_2(0.0);
  auto heads = xit->GetHeads();
  XlinkHead *head0, *head1;
  int ihead = -1;
  bool unbind[2] = {false, false};
  // Do explicitly for now
  head0 = &(*(heads->begin()));
  head1 = &(*(heads->begin()+1));

  // Get all relevant information, need for force dep
  // velocity, etc
  auto aid0 = head0->GetAttach().first;
  auto aidx0 = (*oid_position_map_)[aid0];
  auto rod0 = (*simples_)[aidx0];
  auto crosspos0 = head0->GetAttach().second;
  auto lrod0 = rod0->GetRigidLength();
  auto flink = head0->GetForce();

  auto aid1 = head1->GetAttach().first;
  auto aidx1 = (*oid_position_map_)[aid1];
  auto rod1 = (*simples_)[aidx1];
  auto crosspos1 = head1->GetAttach().second;
  auto lrod1 = rod1->GetRigidLength();

  // Only need one of the forces, should be equal and opposite
  double ui_dot_f = dot_product(ndim_, rod0->GetRigidOrientation(), flink);
  double uj_dot_f = dot_product(ndim_, rod1->GetRigidOrientation(), flink);
  double ui_dot_uj = dot_product(ndim_, rod0->GetRigidOrientation(), rod1->GetRigidOrientation());

  //std::cout << "ui_dot_f: " << std::setprecision(16) << ui_dot_f << std::endl;
  //std::cout << "uj_dot_f: " << std::setprecision(16) << uj_dot_f << std::endl;
  //std::cout << "ui_dot_uj: " << std::setprecision(16) << ui_dot_uj << std::endl;

  double f_mag_i = 0.0, f_mag_j = 0.0;

  if (stall_type_) {
    f_mag_i = ui_dot_f;
    f_mag_j = -uj_dot_f;
  }

  //std::cout << "initial fmagi: " << std::setprecision(16) << f_mag_i << ", fmagj: " << f_mag_j << std::endl;

  // Calculate the velocity scale based on directionality
  double velocity_i = velocity_[0] *
    ((ui_dot_uj > velocity_switch_costheta_[0]) ?
     velocity_p_scale_[0] :
     velocity_ap_scale_[0]);
  double elong = sqrt(dot_product(ndim_, xit->GetRcross(), xit->GetRcross()));
  if (f_mag_i*velocity_i > 0.0) {
    f_mag_i = 0.0;
  } else {
    if (stall_type_ == 1) {
      f_mag_i = ABS(f_mag_i);
    } else if (stall_type_ == 2) {
      f_mag_i = k_stretch_ * elong;
    } else if (stall_type_ == 3) {
      f_mag_i = -ABS(f_mag_i);
    }
  }

  double velocity_j = velocity_[1] *
    ((ui_dot_uj > velocity_switch_costheta_[1]) ?
     velocity_p_scale_[1] :
     velocity_ap_scale_[1]);
  if (f_mag_j*velocity_j > 0.0) {
    f_mag_j = 0.0;
  } else {
    if (stall_type_ == 1)
      f_mag_j = ABS(f_mag_j);
    else if (stall_type_ == 2)
      f_mag_j = k_stretch_ * elong;
    else if (stall_type_ == 3) 
      f_mag_j = -ABS(f_mag_j);
  }

  // Stall force
  if (f_mag_i < f_stall_[0]) {
    if (debug_trace) {
      std::cout << std::setprecision(16) << "Stalling velocity_i: " << velocity_i << " -> ";
    }
    velocity_i *= 1.0 - f_mag_i/f_stall_[0];
    if (debug_trace) {
      std::cout << std::setprecision(16) << velocity_i << std::endl;
    }
  } else {
    velocity_i = 0.0;
  }

  if (f_mag_j < f_stall_[1]) {
    if (debug_trace) {
      std::cout << std::setprecision(16) << "Stalling velocity_j: " << velocity_j << " -> ";
    }
    velocity_j *= 1.0 - f_mag_j/f_stall_[1];
    if (debug_trace) {
      std::cout << std::setprecision(16) << velocity_j << std::endl;
    }
  } else {
    velocity_j = 0.0;
  }

  // Head 0
  head0->SetNExp_1_2(0.0);

  // Update velocity
  crosspos0 += velocity_i * head0->GetDelta();
  if (diffusion_bound_2_[0] > 0.0) {
    crosspos0 += sqrt(2.0 * head0->GetDelta() * diffusion_bound_2_[0]) *
        gsl_ran_gaussian_ziggurat(head0->GetRNG()->r, 1.0);
    crosspos0 += ui_dot_f * diffusion_bound_2_[0] * head0->GetDelta();
  }
  if (crosspos0 > lrod0) {
    crosspos0 = lrod0;
    if (!end_pause_[0]) {
      unbind[0] = true;
    }
  } else if (crosspos0 < 0.0) {
    crosspos0 = 0.0;
    if (!end_pause_[0]) {
      unbind[0] = true;
    }
  }
  head0->Attach(aid0, crosspos0);

  // Head 1
  head1->SetNExp_1_2(0.0);

  // Update velocity
  crosspos1 += velocity_j * head1->GetDelta();
  if (diffusion_bound_2_[1] > 0.0) {
    crosspos1 += sqrt(2.0 * head1->GetDelta() * diffusion_bound_2_[1]) *
        gsl_ran_gaussian_ziggurat(head1->GetRNG()->r, 1.0);
    crosspos1 += -uj_dot_f * diffusion_bound_2_[1] * head1->GetDelta();
  }
  if (crosspos1 > lrod1) {
    crosspos1 = lrod1;
    if (!end_pause_[1]) {
      unbind[1] = true;
    }
  } else if (crosspos1 < 0.0) {
    crosspos1 = 0.0;
    if (!end_pause_[1]) {
      unbind[1] = true;
    }
  }
  head1->Attach(aid1, crosspos1);

  // Run the position update
  xit->UpdateStagePosition(rod0->GetRigidPosition(),
                           rod0->GetRigidOrientation(),
                           rod0->GetRigidLength(),
                           rod0->GetOID(),
                           rod1->GetRigidPosition(),
                           rod1->GetRigidOrientation(),
                           rod1->GetRigidLength(),
                           rod1->GetOID());



  // Check for detachment based on unbind
  if (unbind[0] && unbind[1]) {
    // Double detach
    Detach_2_0(xit);
  } else if (unbind[0] || unbind[1]) {
    // Single detachment
    int whichhead = unbind[0] ? 0 : 1;
    Detach_2_1(xit, whichhead);
  } else {
    //xit->SetPosition(avgpos);
    //xit->SetPrevPosition(oldxitpos);
    //xit->UpdatePeriodic();
    // Figure out which ehads are attached
    nbound2_[0]++;
    nbound2_[1]++;
  }
}

void XlinkKMC::TransferForces() {
  // Transfer the forces from doubly bound xlinks to their
  // bound rods
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Grab the species and start doing the transfer
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    if ((*xit)->GetBoundState() == doubly) {
      ApplyStage2Force(*xit);
    }
  }
}

void XlinkKMC::ApplyStage2Force(Xlink *xit) {
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;

  auto aid0 = head0->GetAttach().first;
  auto aid1 = head1->GetAttach().first;
  auto aidx0 = (*oid_position_map_)[aid0];
  auto aidx1 = (*oid_position_map_)[aid1];
  auto mrod0 = (*simples_)[aidx0];
  auto mrod1 = (*simples_)[aidx1];
  auto rx0 = head0->GetRigidPosition();
  auto rx1 = head1->GetRigidPosition();
  auto sx0 = head0->GetScaledPosition();
  auto sx1 = head1->GetScaledPosition();

  if (debug_trace) {
    std::cout << "[" << head0->GetOID() << ":" << mrod0->GetOID() << "] <-> [" << head1->GetOID() << ":" << mrod1->GetOID()
      << "] Transferring 2stage force\n";
  }

  auto crosspos0 = head0->GetAttach().second;
  auto crosspos1 = head1->GetAttach().second;
  auto lrod0 = mrod0->GetRigidLength();
  auto lrod1 = mrod1->GetRigidLength();
  auto urod0 = mrod0->GetRigidOrientation();
  auto urod1 = mrod1->GetRigidOrientation();

  // Make sure that forces are equal and opposite
  // flink0 is the main one, standing in the place of flink
  double zerovec[3] = {0.0, 0.0, 0.0};
  double flink0[3] = {0.0, 0.0, 0.0};
  double flink1[3] = {0.0, 0.0, 0.0};
  std::copy(head0->GetForce(), head0->GetForce()+ndim_, flink0);
  std::copy(head1->GetForce(), head1->GetForce()+ndim_, flink1);

  double lambda = crosspos0 - 0.5 * lrod0;
  double mu     = crosspos1 - 0.5 * lrod1;

  double rcontact_i[3] = {0.0, 0.0, 0.0};
  double rcontact_j[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < ndim_; ++i) {
    rcontact_i[i] = urod0[i] * lambda;
    rcontact_j[i] = urod1[i] * mu;
  }

  double tau[3];
  double taubond0[3] = {0.0, 0.0, 0.0};
  double taubond1[3] = {0.0, 0.0, 0.0};
  cross_product(rcontact_i, flink0, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
    taubond0[i] += tau[i];
  }
  cross_product(rcontact_j, flink0, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
    taubond1[i] -= tau[i];
  }

  // Apply the force
  mrod0->AddForce(flink0);
  mrod0->AddTorque(taubond0);
  mrod1->AddForce(flink1);
  mrod1->AddTorque(taubond1);

  // Zero the force on the xlink
  head0->SetForce(zerovec);
  head1->SetForce(zerovec);

}

void XlinkKMC::Dump() {
  // print out the information appropriate to kmc
  if (debug_trace) {
    XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
    printf("XlinkKMC -> dump\n");
    printf("\t{n_exp_0_1: %2.4f, n_exp_1_2: %2.4f}\n", pxspec->GetNExp_0_1(), pxspec->GetNExp_1_2());
    printf("\t{nfree:  %d}\n", pxspec->GetNFree());
    printf("\t{nbound1: %d,%d}\n", pxspec->GetNBound1()[0], pxspec->GetNBound1()[1]);
    printf("\t{nbound2: %d,%d}\n", pxspec->GetNBound2()[0], pxspec->GetNBound2()[1]);
    pxspec->DumpKMC();
  }
}

void XlinkKMC::PrepOutputs() {
  kmc_file_name_ << "sc-kmc-XlinkKMC.log";
  kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out);
  kmc_file << "step #ntot #nfree #nbound1[0] #nbound1[1] #nbound2[0] #nbound2[1] nexp01 nexp12\n";
  kmc_file.close();
}

void XlinkKMC::WriteEvent(const std::string &pString) {
  if (write_event_) {
    kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out | std::ios_base::app);
    kmc_file << pString << std::endl;
    kmc_file.close();
  }
}

void XlinkKMC::WriteOutputs(int istep) {
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out | std::ios_base::app);
  kmc_file.precision(16);
  kmc_file.setf(std::ios::fixed, std::ios::floatfield);
  kmc_file << istep << " ";
  kmc_file << pxspec->GetNMembers() << " " << pxspec->GetNFree() << " " <<
    pxspec->GetNBound1()[0] << " " << pxspec->GetNBound1()[1] << " " << 
    pxspec->GetNBound2()[0] << " " << pxspec->GetNBound2()[1] << " ";
  kmc_file << pxspec->GetNExp_0_1() << " "  << pxspec->GetNExp_1_2() << "\n";
  kmc_file.close();
}
