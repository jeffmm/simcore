#include "argon.h"

void Argon::Init() {
  Simple::Init();
  for (int i=0; i<n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = 20*(gsl_rng_uniform_pos(rng_.r)-0.5);
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  orientation_[0] = 0; // Set default color to red
}

void Argon::InitConfigurator(std::array<double, 3> rx, std::array<double, 3> vx) {
  Simple::Init();
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] = 1.0/sqrt(n_dim_);
    velocity_[i] = vx[i];
    position_[i] = rx[i];
    prev_position_[i] = position_[i] - delta_ * velocity_[i];
  }
  UpdatePeriodic();
}

void Argon::UpdatePosition() {
  ZeroForce();
  ApplyInteractions();
  Integrate();
  UpdatePeriodic();
  ClearInteractions();
}

// Update the position based on MP
void Argon::UpdatePositionMP() {
    Integrate();
    UpdatePeriodic();
}

// Basic verlet integrator, very stable
void Argon::Integrate() {
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

void Argon::UpdateKineticEnergy() {
  double vel_mag_sqr = 0.0;
  for (int i=0; i<n_dim_; ++i)
    vel_mag_sqr += SQR(velocity_[i]);
  k_energy_ = 0.5 * mass_ * vel_mag_sqr;
}

double const Argon::GetKineticEnergy() {
  UpdateKineticEnergy();
  return k_energy_;
}

// Species Specific
void ArgonSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "Argon species\n";

  YAML::Node node = YAML::LoadFile(filename);

  // See what kind of insertion we are doing
  std::string insertion_type;
  insertion_type = node["argon"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["argon"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap: " << (can_overlap ? "true" : "false") << std::endl;

  if (insertion_type.compare("xyz") == 0) {
    std::cout << "Nope, not yet!\n";
    exit(1);
  } else if (insertion_type.compare("dat") == 0) {
    int nargon      = node["argon"]["ar"]["num"].as<int>(); 
    double mass     = node["argon"]["ar"]["mass"].as<double>();
    double diameter = node["argon"]["ar"]["diameter"].as<double>();
    std::cout << "   n_argons: " << nargon << std::endl;
    std::cout << "   mass:     " << mass << std::endl;
    std::cout << "   diameter: " << diameter << std::endl;

    // Grab the restart file
    std::string restfile = node["argon"]["ar"]["restfile"].as<std::string>();
    FILE *rfile = fopen(restfile.c_str(), "r");

    std::vector<std::array<double, 3>> positions;
    std::vector<std::array<double, 3>> velocities;
    for (int i = 0; i < nargon; ++i) {
      std::array<double, 3> rx = {0};
      fscanf(rfile, "%lf%lf%lf", &rx[0], &rx[1], &rx[2]);
      positions.push_back(rx);
    }
    for (int i = 0; i < nargon; ++i) {
      std::array<double, 3> vx = {0};
      fscanf(rfile, "%lf%lf%lf", &vx[0], &vx[1], &vx[2]);
      velocities.push_back(vx);
    }
    fclose(rfile);

    for (int i = 0; i < nargon; ++i) {
      Argon *member = new Argon(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(positions[i], velocities[i]); 
      members_.push_back(member);
    }


  } else if (insertion_type.compare("random") == 0) {
    int nmdbeads    = node["argon"]["mdbead"]["num"].as<int>();
    double diameter = node["argon"]["mdbead"]["diameter"].as<double>();
    double mass     = node["argon"]["mdbead"]["mass"].as<double>();
    std::cout << "   n_argons: " << nmdbeads << std::endl;
    std::cout << "   diameter:   " << diameter << std::endl;
    std::cout << "   mass:       " << mass << std::endl;
    params_->n_argon = nmdbeads;
    params_->argon_diameter = diameter;
    params_->argon_mass = mass;

    for (int i = 0; i < nmdbeads; ++i) {
      Argon *member = new Argon(params_, space_, gsl_rng_get(rng_.r), GetSID());
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
            interactionmindist idm;
            auto part1 = member;
            auto part2 = (*mdit);
            MinimumDistance(part1, part2, idm, space_->n_dim, space_->n_periodic, space_);
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
