#include "dynamic_instability.h"


void DynamicInstabilityKMC::Init(space_struct *pSpace,
                                 ParticleTracking *pTracking,
                                 PotentialManager *pPotentials,
                                 SpeciesBase *spec1,
                                 SpeciesBase *spec2,
                                 int ikmc,
                                 YAML::Node &node,
                                 long seed) {
  KMCBase::Init(pSpace, pTracking, pPotentials, spec1, spec2, ikmc, node, seed);

  delta_ = spec1_->GetDelta();

  // Grab our specific claims
  f_shrink_to_grow_   = node["kmc"][ikmc]["f_shrink_to_grow"].as<double>();
  f_shrink_to_pause_  = node["kmc"][ikmc]["f_shrink_to_pause"].as<double>();
  f_pause_to_shrink_  = node["kmc"][ikmc]["f_pause_to_shrink"].as<double>();
  f_pause_to_grow_    = node["kmc"][ikmc]["f_pause_to_grow"].as<double>();
  f_grow_to_shrink_   = node["kmc"][ikmc]["f_grow_to_shrink"].as<double>();
  f_grow_to_pause_    = node["kmc"][ikmc]["f_grow_to_pause"].as<double>();
  // Convert to probabilities
  p_stg_ = f_shrink_to_grow_ * delta_;
  p_stp_ = f_shrink_to_pause_ * delta_;
  p_gts_ = f_grow_to_shrink_ * delta_;
  p_gtp_ = f_grow_to_pause_ * delta_;
  p_ptg_ = f_pause_to_grow_ * delta_;
  p_pts_ = f_pause_to_shrink_ * delta_;
  v_poly_             = node["kmc"][ikmc]["v_poly"].as<double>();
  v_depoly_           = node["kmc"][ikmc]["v_depoly"].as<double>();

  BrRodSpecies* prspec = dynamic_cast<BrRodSpecies*>(spec1_);
  max_length_ = prspec->GetMaxLength();
  min_length_ = prspec->GetMinLength();
}

void DynamicInstabilityKMC::Print() {
  std::cout << "BrRod Dynamic Instability KMC Module\n";
  KMCBase::Print();
  std::cout << std::setprecision(16) << "\tdelta:             " << delta_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_shrink_to_grow:  " << f_shrink_to_grow_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_shrink_to_pause: " << f_shrink_to_pause_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_pause_to_grow:   " << f_pause_to_grow_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_pause_to_shrink: " << f_pause_to_shrink_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_grow_to_pause:   " << f_grow_to_pause_ << std::endl;
  std::cout << std::setprecision(16) << "\tf_grow_to_shrink:  " << f_grow_to_shrink_ << std::endl;
  std::cout << std::setprecision(16) << "\tv_poly:            " << v_poly_ << std::endl;
  std::cout << std::setprecision(16) << "\tv_depoly:          " << v_depoly_ << std::endl;
  std::cout << std::setprecision(16) << "\tmax_length:        " << max_length_ << std::endl;
  std::cout << std::setprecision(16) << "\tmin_length:        " << min_length_ << std::endl;
}

void DynamicInstabilityKMC::Dump() {
  if (debug_trace) {
    std::cout << "DynamicInstabilityKMC -> dump\n";
  }
}

void DynamicInstabilityKMC::StepKMC() {
  UpdatePolymerizationState();
  GrowBonds();
}

void DynamicInstabilityKMC::UpdatePolymerizationState() {
  BrRodSpecies *prspec = dynamic_cast<BrRodSpecies*>(spec1_);
  //std::cout << "UpdatePolymerizationState\n";
  auto rod_members = prspec->GetMembers();

  for (auto rodit = rod_members->begin(); rodit != rod_members->end(); ++rodit) {
    BrRod *rod = *rodit;
    //std::cout << "Examining rod: " << rod->GetOID() << ", rid: " << rod->GetRID() << ", cid: " << rod->GetCID() << std::endl;
    //rod->Dump();

    double f_stabilize_fr = 1.0;
    double f_stabilize_fc = 1.0;
    double f_stabilize_vg = 1.0;
    double f_stabilize_vs = 1.0;
    int stab_state = rod->GetStabilizationState(&f_stabilize_fr, &f_stabilize_fc, &f_stabilize_vg, &f_stabilize_vs);
    //std::cout << "stab state: " << stab_state << std::endl;
    //std::cout << "f_stabilize_fr: " << std::setprecision(16) << f_stabilize_fr << std::endl;
    //std::cout << "f_stabilize_fc: " << std::setprecision(16) << f_stabilize_fc << std::endl;
    poly_state_t poly_state = rod->GetPolyState();
    poly_state_t poly_state_new = poly_state;
    //std::cout << "poly state: " << poly_state << std::endl;
    auto mrng = rod->GetRNG();

    // Adjust DI rates if MT is stabilized
    double p_factor_r = 1.0;
    double p_factor_c = 1.0;

    if (stab_state == 1) {
      p_factor_r = f_stabilize_fr;
      p_factor_c = f_stabilize_fc;
    }

    double p_stg = p_stg_ * p_factor_r;
    double p_stp = p_stp_ * p_factor_r;
    double p_gts = p_gts_ * p_factor_c;
    double p_gtp = p_gtp_ * p_factor_c;
    double p_ptg = p_ptg_ * p_factor_r;
    double p_pts = p_pts_ * p_factor_c;

    //std::cout << std::setprecision(16) << "p_stg: " << p_stg << std::endl
    //  << "p_stp: " << p_stp << std::endl
    //  << "p_gts: " << p_gts << std::endl
    //  << "p_gtp: " << p_gtp << std::endl
    //  << "p_ptg: " << p_ptg << std::endl
    //  << "p_pts: " << p_pts << std::endl;

    double p_norm = 0.0;

    // XXX FIXME Check the wall potential for the force induced catstrophe stuff here

    // Can now start checking states
    // Shrinking state
    double roll = gsl_rng_uniform_pos(mrng->r);
    //std::cout << "roll: " << std::setprecision(16) << roll << std::endl;
    if (poly_state == SHRINK) {
      p_norm = p_stg + p_stp;
      //std::cout << "Shrink p_norm: " << std::setprecision(16) << p_norm << std::endl;
      if (p_norm > 1.0) {
        if (roll < p_stg/p_norm) {
          poly_state_new = GROW;
        } else {
          poly_state_new = PAUSE;
        }
      } else {
        if (roll < p_stg) {
          poly_state_new = GROW;
        } else if (roll >= p_stg && roll < (p_stg + p_stp)) {
          poly_state_new = PAUSE;
        }
      }
    } // shrinking

    // MT growing
    else if (poly_state == GROW) {
      p_norm = p_gts + p_gtp;
      //std::cout << "Grow p_norm: " << std::setprecision(16) << p_norm << std::endl;
      if (p_norm > 1.0) {
        if (roll < p_gts/p_norm) {
          poly_state_new = SHRINK;
        } else {
          poly_state_new = PAUSE;
        }
      } else {
        if (roll < p_gts) {
          poly_state_new = SHRINK;
        } else if (roll >= p_gts && roll < (p_gts + p_gtp)) {
          poly_state_new = PAUSE;
        }
      }
    } // growing

    // MT paused
    else if (poly_state == PAUSE) {
      p_norm = p_ptg + p_pts;
      //std::cout << "Pause p_norm: " << std::setprecision(16) << p_norm << std::endl;
      if (p_norm > 1.0) {
        if (roll < p_ptg/p_norm) {
          poly_state_new = GROW;
        } else {
          poly_state_new = SHRINK;
        }
      } else {
        if (roll < p_ptg) {
          poly_state_new = GROW;
        } else if (roll >= p_ptg && roll < (p_ptg + p_pts)) {
          poly_state_new = PAUSE;
        }
      }
    } // end of MT poly check

    // Final check if we hit max or min lengths
    if (rod->GetLength() > max_length_) {
      poly_state_new = SHRINK;
    } else if (rod->GetLength() <= min_length_) {
      poly_state_new = GROW;
    }

    if (poly_state != poly_state_new) {
      rod->SetPolyState(poly_state_new);
      if (debug_trace) {
        // XXX FIXME write this event
        std::cout << "Changing from state: " << poly_state << " to " << poly_state_new << std::endl;
      }
    }
  } // loop over rods

}

void DynamicInstabilityKMC::GrowBonds() {
  BrRodSpecies *prspec = dynamic_cast<BrRodSpecies*>(spec1_);
  //std::cout << "GrowBonds\n";
  auto rod_members = prspec->GetMembers();

  double delta_L_depoly = -v_depoly_ * delta_;
  double delta_L_poly = v_poly_ * delta_;

  for (auto rodit = rod_members->begin(); rodit != rod_members->end(); ++rodit) {
    BrRod *rod = *rodit;
    //std::cout << "Grow rod: " << rod->GetOID() << ", rid: " << rod->GetRID() << ", cid: " << rod->GetCID() << std::endl;
    //rod->Dump();

    poly_state_t poly_state = rod->GetPolyState();
    if (poly_state == PAUSE) continue; // dont' do anything on a pause

    // XXX FIXME Check the wall potentials for vp dep and stuff

    double delta_L = 0.0;
    if (poly_state == GROW)
      delta_L = delta_L_poly;
    else if (poly_state == SHRINK)
      delta_L = delta_L_depoly;
    else if (poly_state == PAUSE)
      delta_L = 0.0;
    else
      delta_L = 0.0;

    //std::cout << "Delta_L: " << std::setprecision(16) << delta_L << std::endl;
    // Check stabilizations
    double f_stabilize_fr = 1.0;
    double f_stabilize_fc = 1.0;
    double f_stabilize_vg = 1.0;
    double f_stabilize_vs = 1.0;
    auto stab_state = rod->GetStabilizationState(&f_stabilize_fr, &f_stabilize_fc, &f_stabilize_vg, &f_stabilize_vs);
    //std::cout << "stab state: " << poly_state << std::endl;
    //std::cout << "f_stabilize_vg: " << std::setprecision(16) << f_stabilize_vg << std::endl;
    //std::cout << "f_stabilize_vs: " << std::setprecision(16) << f_stabilize_vs << std::endl;
    if (stab_state == 1) {
      // Modify the growth and shrink velocities
      if (poly_state == GROW) {
        delta_L = delta_L * f_stabilize_vg;
      } else if (poly_state == PAUSE) {
        delta_L = delta_L * f_stabilize_vs;
      }
    }

    // Push this info to the rod so it can update child bonds, etc
    rod->UpdateRodLength(delta_L);
  }
}

