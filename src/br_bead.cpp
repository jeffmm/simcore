#include "br_bead.h"

void BrBead::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}

void BrBead::UpdatePosition() {
  KickBead();
  ApplyInteractions();
  for (int i=0; i<n_dim_; ++i)
    position_[i] = position_[i] + force_[i] * delta_ / diameter_;
  UpdatePeriodic();
  ClearInteractions();
  ZeroForce();
}

void BrBead::UpdatePositionMP() {
    KickBead();
    std::copy(position_, position_+n_dim_, prev_position_);
    for (int i = 0; i < n_dim_; ++i) {
        position_[i] = position_[i] + force_[i] * delta_ / diameter_;
        dr_tot_[i] += position_[i] - prev_position_[i];
    }
    UpdatePeriodic();
}

//void BrBeadSpecies::InitPotentials (system_parameters *params) {
  //AddPotential(SID::br_bead, SID::br_bead, 
      //new WCA(params->lj_epsilon,params->br_bead_diameter,
                //space_, pow(2.0, 1.0/6.0)*params->br_bead_diameter));
//}

void BrBeadSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "BRBead species\n";

  YAML::Node node = YAML::LoadFile(filename);

  // See what kind of insertion we are doing
  std::string insertion_type;
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  insertion_type = node["br_bead"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["br_bead"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap: " << (can_overlap ? "true" : "false") << std::endl;
  if (node["br_bead"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["br_bead"]["properties"]["color"][i].as<double>();
    }
  }
  std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
    << color[3] << "]\n";

  if (insertion_type.compare("xyz") == 0) {
    std::cout << "Nope, not yet!\n";
    exit(1);
  } else if (insertion_type.compare("random") == 0) {
    int nbrbeads    = node["br_bead"]["brbead"]["num"].as<int>();
    double diameter = node["br_bead"]["brbead"]["diameter"].as<double>();
    std::cout << "   n_br_beads: " << nbrbeads << std::endl;
    std::cout << "   diameter:   " << diameter << std::endl;
    params_->n_br_bead = nbrbeads;
    params_->br_bead_diameter = diameter;

    for (int i = 0; i < nbrbeads; ++i) {
      BrBead *member = new BrBead(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      member->SetColor(color, 0);
      
      // Check if overlaps allowed
      if (can_overlap) {
        members_.push_back(member);
      } else {
        // Check against other br beads
        bool isoverlap = true;
        int numoverlaps = 0;
        do {
          numoverlaps++;
          isoverlap = false;
          for (auto brit = members_.begin(); brit != members_.end() && !isoverlap; ++brit) {
            interactionmindist idm;
            auto part1 = member;
            auto part2 = (*brit);
            MinimumDistance(part1, part2, idm, space_->n_dim, space_->n_periodic, space_);
            double diameter2 = diameter*diameter;

            if (idm.dr_mag2 < diameter2) {
              isoverlap = true;
              member->Init();
            }
          }
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << nbrbeads;
            std::cout << ".  Check packing ratio for objects.\n";
            exit(1);
          }
        } while (isoverlap);
        members_.push_back(member);
      }
    }

  } else {
    std::cout << "Nope, not yet!\n";
    exit(1);
  }
}



