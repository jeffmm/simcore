#include "dy_rod.h"
#include <yaml-cpp/yaml.h>

void DyRod::Init() {
  if (diffusion_validation_flag_) {
    for (int i=0; i<n_dim_; ++i)
      position_[i] = orientation_[i] = 0.0;
    orientation_[n_dim_-1] = 1.0;
    UpdatePeriodic();
  }
  else {
    InsertRandom(0.5*length_+diameter_);
  }
  SetDiffusion();
  std::fill(body_frame_, body_frame_+6, 0.0);
  // Init bond lengths and diameter
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetRigidPosition(position_);
    bond->SetRigidScaledPosition(scaled_position_);
    bond->SetRigidOrientation(orientation_);
    bond->SetRigidLength(length_);
    bond->SetRigidDiameter(diameter_);
    bond->SetLength(child_length_);
    bond->SetDiameter(diameter_);
  }
  // Set positions for sites and bonds
  UpdatePeriodic();
  UpdateBondPositions();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)
    bond->UpdatePeriodic();
}

void DyRod::InitConfigurator(const double* const x, const double* const u, const double l) {
  length_ = l;
  SetPosition(x);
  SetOrientation(u);
  UpdatePeriodic();
  SetDiffusion();
  std::fill(body_frame_, body_frame_+6, 0.0);
  // Init bond lengths and diameter
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetRigidPosition(position_);
    bond->SetRigidScaledPosition(scaled_position_);
    bond->SetRigidOrientation(orientation_);
    bond->SetRigidLength(length_);
    bond->SetRigidDiameter(diameter_);
    bond->SetLength(child_length_);
    bond->SetDiameter(diameter_);
  }
  // Set positions for sites and bonds
  UpdatePeriodic();
  UpdateBondPositions();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)
    bond->UpdatePeriodic();
  UpdateSitePositions();
}

void DyRod::ApplyForcesTorques() {
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    AddForce(bond->GetForce());
    AddTorque(bond->GetTorque());
    AddPotential(bond->GetPotentialEnergy());
  }
  // Add driving
  double f_dr[3];
  for (int i=0; i<n_dim_; ++i)
    f_dr[i] = orientation_[i]*driving_factor_;
  AddForce(f_dr);
}

void DyRod::UpdatePositionMP() {
  ApplyForcesTorques();
  Integrate();
  UpdateSitePositions();
  // Update end site positions for tracking trajectory for neighbors
  UpdateBondPositions();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) 
    bond->UpdatePeriodic();
}

void DyRod::UpdateSitePositions() {
  // First set prev positions for sites
  elements_[0].SetPrevPosition(elements_[0].GetPosition());
  elements_[1].SetPrevPosition(elements_[1].GetPosition());
  // Then update site positions
  double pos[3];
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] - 0.5 * length_ * orientation_[i];
  elements_[0].SetPosition(pos);
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] + 0.5 * length_ * orientation_[i];
  elements_[1].SetPosition(pos);
  // update trajectories of end sites for neighbor lists
  elements_[0].AddDr();
  elements_[1].AddDr();
  // Update positions based on PBCs for rod
  UpdatePeriodic();
  // Record new site positions to avoid issues with pbcs for trajectories
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] - 0.5 * length_ * orientation_[i];
  elements_[0].SetPosition(pos);
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] + 0.5 * length_ * orientation_[i];
  elements_[1].SetPosition(pos);
}

void DyRod::UpdateBondPositions() {
  // Set site of first bond COM and update remaining COMs
  double pos[3];
  for (int i=0; i<n_dim_; ++i)
    pos[i] = position_[i] + 0.5*(child_length_-length_)*orientation_[i];
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    bond->SetRigidPosition(position_);
    bond->SetRigidLength(length_);
    bond->SetRigidDiameter(diameter_);
    bond->SetRigidScaledPosition(scaled_position_);
    bond->SetRigidOrientation(orientation_);
    bond->SetPosition(pos);
    bond->SetOrientation(orientation_);
    // Set next bond COM
    for (int i=0; i<n_dim_; ++i)
      pos[i] += orientation_[i] * child_length_;
  }
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
void DyRod::Integrate() {
  //Explicit calculation of Xi.F_s
  for (int i=0; i<n_dim_; ++i) {
    for (int j=0; j<n_dim_; ++j) {
      position_[i] +=
        gamma_par_*orientation_[i]*orientation_[j]*force_[j]*delta_;
    }
    position_[i] += force_[i]*gamma_perp_*delta_;
  }
  //Add the random displacement dr(t)
  AddRandomDisplacement();
  //Update the orientation due to torques and random rotation
  UpdateOrientation();
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void DyRod::AddRandomDisplacement() {
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_par_);
  for (int i=0; i<n_dim_; ++i)
    position_[i] += mag * orientation_[i];
  // Then the perpendicular component(s)
  for (int j=0; j<n_dim_-1; ++j) {
    mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_perp_);
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
void DyRod::UpdateOrientation() {
  // First handle reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, 3); // ndim=3 since torques
  for (int i=0; i<n_dim_; ++i)
    orientation_[i] += du[i]*delta_/gamma_rot_;
  // Now handle the random orientation update
  // If we aren't supposed to diffuse, just return
    for (int j=0; j<n_dim_-1; ++j) {
      double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
      for (int i=0; i<n_dim_; ++i)
        orientation_[i] += mag * body_frame_[n_dim_*j+i];
    }
  normalize_vector(orientation_, n_dim_);
}

/* calculates vector(s) orthogonal to orientation of rod */
void DyRod::GetBodyFrame() {
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

/* Initialize diffusion coefficients and std dev for random numbers */
void DyRod::SetDiffusion() {
  double logLD = log(length_/diameter_);
  gamma_par_ = 2.0*length_ / (3.0*logLD);
  gamma_perp_ = 2.0*gamma_par_;
  gamma_rot_ = length_*length_*length_ / (9.0*logLD);
  rand_sigma_par_ = sqrt(2.0*delta_/gamma_par_);
  rand_sigma_perp_ = sqrt(2.0*delta_/gamma_perp_);
  rand_sigma_rot_ = sqrt(2.0*delta_/gamma_rot_);
}

void DyRod::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)  {
    bond->SetColor(color_, draw_type_);
    bond->Draw(graph_array);
  }
}

void DyRod::Dump() {
  printf("{%d,%d,%d} -> ", GetOID(), GetRID(), GetCID());
  printf("x(%2.2f, %2.2f, %2.2f), ", GetPosition()[0], GetPosition()[1], GetPosition()[2]);
  printf("u(%2.2f, %2.2f, %2.2f), ", orientation_[0], orientation_[1], orientation_[2]);
  printf("f(%2.2f, %2.2f, %2.2f), ", GetForce()[0], GetForce()[1], GetForce()[2]);
  printf("t(%2.2f, %2.2f, %2.2f), ", GetTorque()[0], GetTorque()[1], GetTorque()[2]);
  printf("l(%2.2f), ke(%2.2f), pe(%2.2f)\n", length_, GetKineticEnergy(), GetPotentialEnergy());
  printf("\tsites ->\n");
  for (auto it=elements_.begin(); it!=elements_.end(); ++it) {
    printf("\t");
    it->Dump();
  }
  printf("\tbonds ->\n");
  for (auto it=v_elements_.begin(); it!=v_elements_.end();++it) {
    printf("\t");
    it->Dump();
  }
}

// Species specifics
void DyRodSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "DyRod species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";
  std::string insertion_type;
  insertion_type = node["dy_rod"]["properties"]["insertion_type"].as<std::string>();

  bool can_overlap = node["dy_rod"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap:        " << (can_overlap ? "true" : "false") << std::endl;

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 1; // default to orientation
  if (node["dy_rod"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["dy_rod"]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }

  if (insertion_type.compare("xyz") == 0) {
    if (!can_overlap) {
      std::cout << "Warning, location insertion overrides overlap\n";
      can_overlap = true;
    }
    max_length_ = node["dy_rod"]["properties"]["max_length"].as<double>();
    std::cout << "   max length:     " << max_length_ << std::endl;
    min_length_ = node["dy_rod"]["properties"]["min_length"].as<double>();
    std::cout << "   min length:     " << min_length_ << std::endl;
    int nrods = (int)node["dy_rod"]["rod"].size();
    std::cout << "   nrods: " << nrods << std::endl;
    params_->n_rod = nrods;
    params_->max_rod_length = max_length_;
    params_->min_rod_length = min_length_;
    for (int irod = 0; irod < nrods; ++irod) {
      double x[3] = {0.0, 0.0, 0.0};
      double u[3] = {0.0, 0.0, 0.0};
      double rlength = 0.0;
      x[0] = node["dy_rod"]["rod"][irod]["x"][0].as<double>();
      x[1] = node["dy_rod"]["rod"][irod]["x"][1].as<double>();
      x[2] = node["dy_rod"]["rod"][irod]["x"][2].as<double>();
      std::cout << "   x(" << x[0] << ", " << x[1] << ", " << x[2] << ")\n";
      u[0] = node["dy_rod"]["rod"][irod]["u"][0].as<double>();
      u[1] = node["dy_rod"]["rod"][irod]["u"][1].as<double>();
      u[2] = node["dy_rod"]["rod"][irod]["u"][2].as<double>();
      normalize_vector(u, space_->n_dim);
      std::cout << "   u(" << u[0] << ", " << u[1] << ", " << u[2] << ")\n";
      rlength = node["dy_rod"]["rod"][irod]["length"].as<double>();
      std::cout << "   length[" << rlength << "]\n";

      DyRod *member = new DyRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(x, u, rlength);
      member->SetColor(color, draw_type);
      //member->SetAnchors(anchors_);
      member->Dump();
      members_.push_back(member);
    }
  } else if (insertion_type.compare("random") == 0) {
    int nrods         = node["dy_rod"]["rod"]["num"].as<int>();
    double rlength    = node["dy_rod"]["rod"]["length"].as<double>();
    double max_length = node["dy_rod"]["rod"]["max_length"].as<double>();
    double min_length = node["dy_rod"]["rod"]["min_length"].as<double>();
    double diameter   = node["dy_rod"]["rod"]["diameter"].as<double>();

    //Set fill type for neumatics and other kinds of forms
    std::cout << std::setw(25) << std::left << "   n rods:" << std::setw(10)
      << std::left << nrods << std::endl;
    std::cout << std::setw(25) << std::left << "   length:" << std::setw(10)
      << std::left << rlength << std::endl;
    std::cout << std::setw(25) << std::left << "   max length:" << std::setw(10)
      << std::left << max_length << std::endl;
    std::cout << std::setw(25) << std::left << "   min length:" << std::setw(10)
      << std::left << min_length << std::endl;
    std::cout << std::setw(25) << std::left << "   diameter:" << std::setw(10)
      << std::left << diameter << std::endl;

    params_->n_rod = nrods;
    params_->rod_length = rlength;
    params_->max_rod_length = max_length;
    params_->min_rod_length = min_length;
    max_length_ = max_length;
    min_length_ = min_length;
    params_->rod_diameter = diameter;
   

    for (int i_mem = 0; i_mem < nrods; ++i_mem) {
      DyRod *member = new DyRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->SetColor(color, draw_type);
      member->Init();

        // Check against all other rods in the sytem
      if (can_overlap) {
        // Done, add to members
        members_.push_back(member);
      } else {
        // Check against all other rods in system
        bool isoverlap = true;
        int numoverlaps = 0;
        do {
          numoverlaps++;
          isoverlap = false;
          for (auto rodit = members_.begin(); rodit != members_.end() && !isoverlap; ++rodit) {
            interactionmindist idm;
            // Just check the 0th element of each
            auto part1 = member->GetSimples()[0];
            auto part2 = (*rodit)->GetSimples()[0];
            MinimumDistance(part1, part2, idm, space_->n_dim, space_->n_periodic, space_);
            double diameter2 = diameter*diameter;

            if (idm.dr_mag2 < diameter2) {
              isoverlap = true;
              member->Init();
            }
          } // check against current members
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i_mem << " of " << nrods;
            std::cout << ".  Check packing ratio for objects.\n";
            exit(1);
          }
        } while (isoverlap);
        members_.push_back(member);
      }
    }
  } 
      
  else {
    printf("nope, not yet\n");
    exit(1);
  }

  n_members_ = members_.size();
}

