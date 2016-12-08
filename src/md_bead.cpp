#include "md_bead.h"

void MDBead::Init() {
  Simple::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 20*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  orientation_[0] = 0; // Set default color to red
}

// Update the position based 
void MDBead::UpdatePosition() {
    Integrate();
    UpdatePeriodic();
}

// Basic verlet integrator, very stable
void MDBead::Integrate() {
  double delta2 = SQR(delta_);
  double temp_pos[3];
  for (int i=0; i<n_dim_; ++i) {
    temp_pos[i] = position_[i];
    position_[i] = 2.0*position_[i] - prev_position_[i] 
      + ((force_[i]/mass_))*delta2;
    velocity_[i] =  (position_[i] - prev_position_[i])/(2.0*delta_);
    prev_position_[i] = temp_pos[i];
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
}

void MDBead::UpdateKineticEnergy() {
  double vel_mag_sqr = 0.0;
  for (int i=0; i<n_dim_; ++i)
    vel_mag_sqr += SQR(velocity_[i]);
  k_energy_ = 0.5 * mass_ * vel_mag_sqr;
}

double const MDBead::GetKineticEnergy() {
  UpdateKineticEnergy();
  return k_energy_;
}

// Species Specific
void MDBeadSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "MDBead species\n";

  YAML::Node node = YAML::LoadFile(filename);

  // See what kind of insertion we are doing
  std::string insertion_type;
  insertion_type = node["md_bead"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["md_bead"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap: " << (can_overlap ? "true" : "false") << std::endl;

  if (insertion_type.compare("xyz") == 0) {
    std::cout << "Nope, not yet!\n";
    exit(1);
  } else if (insertion_type.compare("random") == 0) {
    int nmdbeads    = node["md_bead"]["mdbead"]["num"].as<int>();
    double diameter = node["md_bead"]["mdbead"]["diameter"].as<double>();
    double mass     = node["md_bead"]["mdbead"]["mass"].as<double>();
    std::cout << "   n_md_beads: " << nmdbeads << std::endl;
    std::cout << "   diameter:   " << diameter << std::endl;
    std::cout << "   mass:       " << mass << std::endl;
    params_->n_md_bead = nmdbeads;
    params_->md_bead_diameter = diameter;
    params_->md_bead_mass = mass;

    for (int i = 0; i < nmdbeads; ++i) {
      MDBead *member = new MDBead(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      
      // Check if overlaps allowed
      if (can_overlap) {
        members_.push_back(member);
      } else {
        // Check against other md beads
        bool isoverlap = true;
        int numoverlaps = 0;
        do {
          numoverlaps++;
          isoverlap = false;
          for (auto mdit = members_.begin(); mdit != members_.end() && !isoverlap; ++mdit) {
            Interaction idm;
            auto part1 = member;
            auto part2 = (*mdit);
            MinimumDistance(part1, part2, &idm, space_->n_dim, space_->n_periodic, space_);
            double diameter2 = diameter*diameter;

            if (idm.dr_mag2 < diameter2) {
              isoverlap = true;
              member->Init();
            }
          }
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << nmdbeads;
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
