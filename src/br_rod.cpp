#include "br_rod.h"

#include <yaml-cpp/yaml.h>

#include <iomanip>

#include "sphero_overlap.h"

void BrRod::Init() {
  InsertRandom(0.5*length_+diameter_);
  poly_state_ = GROW;
  stabilization_state_ = 0;
  f_stabilize_fr_ = 1.0;
  f_stabilize_fc_ = 1.0;
  f_stabilize_vg_ = 1.0;
  f_stabilize_vs_ = 1.0;
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
  UpdateRodLength(0.0);
  UpdateBondPositions();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)
    bond->UpdatePeriodic();
}

void BrRod::InitConfigurator(const double* const x, const double* const u, const double l) {
  length_ = l;
  SetPosition(x);
  SetOrientation(u);
  UpdatePeriodic();
  poly_state_ = GROW;
  stabilization_state_ = 0;
  f_stabilize_fr_ = 1.0;
  f_stabilize_fc_ = 1.0;
  f_stabilize_vg_ = 1.0;
  f_stabilize_vs_ = 1.0;
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
  UpdateRodLength(0.0);
  UpdateBondPositions();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)
    bond->UpdatePeriodic();
  UpdateSitePositions();
}

void BrRod::ApplyForcesTorques() {
  //ZeroForce()     BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
  //member->Init();;
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) {
    AddForce(bond->GetForce());
    AddTorque(bond->GetTorque());
    AddPotential(bond->GetPotentialEnergy());
  }
}

void BrRod::UpdatePositionMP() {
  ApplyForcesTorques();
  Integrate();
  //UpdatePeriodic();
  UpdateSitePositions();
  UpdateRodLength(0.0);
  // Update end site positions for tracking trajectory for neighbors
  UpdateBondPositions();
  UpdateAnchors();
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) 
    bond->UpdatePeriodic();
}

void BrRod::UpdateAnchors() {
  // find the anchors
  if (anchors_ == nullptr) return;

  // Update our anchor relative position
  bool found = false;
  for (auto ait = anchors_->begin(); ait != anchors_->end() && !found; ++ait) {
    // Look for our entry in the vector
    for (auto p = ait->second.begin(); p != ait->second.end(); ++p) {
      if (p->idx_other_ == GetSimples()[0]->GetOID()) {

        auto mrod = GetSimples()[0];
        //std::cout << "Old Anchor pos1: (" << std::setprecision(16)
        //  << p->pos1_[0] << ", " << p->pos1_[1] << ", " << p->pos1_[2] << ")\n";
        //std::cout << "Old Anchor pos_rel1: (" << std::setprecision(16)
        //  << p->pos_rel1_[0] << ", " << p->pos_rel1_[1] << ", " << p->pos_rel1_[2] << ")\n";
        for (int i = 0; i < 3; ++i) {
          p->pos1_[i] = mrod->GetRigidPosition()[i] - 0.5 * mrod->GetRigidLength() *
            mrod->GetRigidOrientation()[i];
        }
        p->pos_rel1_[0] = -0.5 * mrod->GetRigidLength();
        p->pos_rel1_[1] = p->pos_rel1_[2] = 0.0;

        //std::cout << "New Anchor pos1: (" << std::setprecision(16)
        //  << p->pos1_[0] << ", " << p->pos1_[1] << ", " << p->pos1_[2] << ")\n";
        //std::cout << "New Anchor pos_rel1: (" << std::setprecision(16)
        //  << p->pos_rel1_[0] << ", " << p->pos_rel1_[1] << ", " << p->pos_rel1_[2] << ")\n";

        found = true;
        break;
      }
    }
  }
}

void BrRod::UpdateSitePositions() {
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

void BrRod::UpdateBondPositions() {
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
void BrRod::Integrate() {
  if (rod_fixed_ == 1) return;
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
void BrRod::AddRandomDisplacement() {
  // Only do this if we are supposed to diffuse
  if (!rod_diffusion_) return;
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
void BrRod::UpdateOrientation() {
  // First handle reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, 3); // ndim=3 since torques
  for (int i=0; i<n_dim_; ++i)
    orientation_[i] += du[i]*delta_/gamma_rot_;
  // Now handle the random orientation update
  // If we aren't supposed to diffuse, just return
  if (rod_diffusion_) {
    for (int j=0; j<n_dim_-1; ++j) {
      double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
      for (int i=0; i<n_dim_; ++i)
        orientation_[i] += mag * body_frame_[n_dim_*j+i];
    }
  }
  normalize_vector(orientation_, n_dim_);
}

/* calculates vector(s) orthogonal to orientation of rod */
void BrRod::GetBodyFrame() {
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
void BrRod::SetDiffusion() {
  double logLD = log(length_/diameter_);
  gamma_par_ = 2.0*length_ / (3.0*logLD);
  gamma_perp_ = 2.0*gamma_par_;
  gamma_rot_ = length_*length_*length_ / (9.0*logLD);
  rand_sigma_par_ = sqrt(2.0*delta_/gamma_par_);
  rand_sigma_perp_ = sqrt(2.0*delta_/gamma_perp_);
  rand_sigma_rot_ = sqrt(2.0*delta_/gamma_rot_);
}

double BrRod::UpdateRodLength(double delta_length) {
  // Intake this information from KMC
  // delta_length is already done via depoly or poly, 
  // so don't have to worry about that here
  length_ += delta_length;
  // Update the bond lengths
  child_length_ = length_/n_bonds_;
  // If necessary, add or remove a bond
  if (child_length_ > max_child_length_) {
    if (debug_trace) {
      std::cout << "UpdateRodLength: Increasing number of bonds from " << n_bonds_
        << " to " << n_bonds_+1 << std::endl;
    }
    n_bonds_++;
    child_length_ = length_/n_bonds_;
    Bond b(v_elements_[0]);
    // Give the new bond a unique OID
    b.InitOID();
    v_elements_.push_back(b);
  } else if (child_length_ < min_length_ && v_elements_.size() > 1) {
    if (debug_trace) {
      std::cout << "UpdateRodLength: Decreasing number of bonds from " << n_bonds_
        << " to " << n_bonds_-1 << std::endl;
    }
    n_bonds_--;
    child_length_ = length_/n_bonds_;
    // We have to get the force and torque, and transfer them onto a different
    // bond
    Bond *firstbond = &(*v_elements_.begin());
    Bond *lastbond = &(*v_elements_.rbegin());
    //std::cout << "firstbond: " << firstbond << ", oid: " << firstbond->GetOID() << std::endl;
    //std::cout << "lastbond: " << lastbond << ", oid: " << lastbond->GetOID() << std::endl;
    firstbond->AddForce(lastbond->GetForce());
    firstbond->AddTorque(lastbond->GetTorque());
    firstbond->AddPotential(lastbond->GetPotentialEnergy());
    v_elements_.pop_back();
  }
  for (auto bond = v_elements_.begin(); bond != v_elements_.end(); ++bond) {
    bond->SetLength(child_length_);
  }

  // Update our diffusion stuff
  SetDiffusion();
  return length_;
}

void BrRod::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond)  {
    bond->SetColor(color_, draw_type_);
    bond->Draw(graph_array);
  }
}

void BrRod::Dump() {
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
void BrRodSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "BrRod species\n";

  YAML::Node node = YAML::LoadFile(filename);

  std::cout << " Generic Properties:\n";

  // See what kind of insertion type, as that changes how we read info
  std::string insertion_type;
  insertion_type = node["br_rod"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["br_rod"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap:        " << (can_overlap ? "true" : "false") << std::endl;

  // Coloring
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  int draw_type = 1; // default to orientation
  if (node["br_rod"]["properties"]["color"]) {
    for (int i = 0; i < 4; ++i) {
      color[i] = node["br_rod"]["properties"]["color"][i].as<double>();
    }
    std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
      << color[3] << "]\n";
  }
  if (node["br_rod"]["properties"]["draw_type"]) {
    std::string draw_type_s = node["br_rod"]["properties"]["draw_type"].as<std::string>();
    std::cout << "   draw_type: " << draw_type_s << std::endl;
    if (draw_type_s.compare("flat") == 0) {
      draw_type = 0;
    } else if (draw_type_s.compare("orientation") == 0) {
      draw_type = 1;
    }
  }

  if (insertion_type.compare("xyz") == 0) {
    if (!can_overlap) {
      std::cout << "Warning, location insertion overrides overlap\n";
      can_overlap = true;
    }
    max_length_ = node["br_rod"]["properties"]["max_length"].as<double>();
    std::cout << "   max length:     " << max_length_ << std::endl;
    min_length_ = node["br_rod"]["properties"]["min_length"].as<double>();
    std::cout << "   min length:     " << min_length_ << std::endl;
    int nrods = (int)node["br_rod"]["rod"].size();
    std::cout << "   nrods: " << nrods << std::endl;
    params_->n_rod = nrods;
    params_->max_rod_length = max_length_;
    params_->min_rod_length = min_length_;
    for (int irod = 0; irod < nrods; ++irod) {
      double x[3] = {0.0, 0.0, 0.0};
      double u[3] = {0.0, 0.0, 0.0};
      double rlength = 0.0;
      x[0] = node["br_rod"]["rod"][irod]["x"][0].as<double>();
      x[1] = node["br_rod"]["rod"][irod]["x"][1].as<double>();
      x[2] = node["br_rod"]["rod"][irod]["x"][2].as<double>();
      std::cout << "   x(" << x[0] << ", " << x[1] << ", " << x[2] << ")\n";
      u[0] = node["br_rod"]["rod"][irod]["u"][0].as<double>();
      u[1] = node["br_rod"]["rod"][irod]["u"][1].as<double>();
      u[2] = node["br_rod"]["rod"][irod]["u"][2].as<double>();
      normalize_vector(u, space_->n_dim);
      std::cout << "   u(" << u[0] << ", " << u[1] << ", " << u[2] << ")\n";
      rlength = node["br_rod"]["rod"][irod]["length"].as<double>();
      std::cout << "   length[" << rlength << "]\n";

      BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->InitConfigurator(x, u, rlength);
      member->SetColor(color, draw_type);
      member->Dump();
      members_.push_back(member);
    }
  } else if (insertion_type.compare("random") == 0) {
    int nrods         = node["br_rod"]["rod"]["num"].as<int>();
    double rlength    = node["br_rod"]["rod"]["length"].as<double>();
    double max_length = node["br_rod"]["rod"]["max_length"].as<double>();
    double min_length = node["br_rod"]["rod"]["min_length"].as<double>();
    double diameter   = node["br_rod"]["rod"]["diameter"].as<double>();

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

    for (int i = 0; i < nrods; ++i) {
      BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      member->SetColor(color, draw_type);

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
              /*if (debug_trace) {
                printf("Overlap detected [oid: %d,%d], [cid: %d, %d] -> (%2.2f < %2.2f)\n",
                        part1->GetOID(), part2->GetOID(), part1->GetCID(), part2->GetCID(), idm.dr_mag2, diameter2);
              }*/
              // We can just call init again to get new random numbers
              member->Init();
            }
          } // check against current members
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << nrods;
            std::cout << ".  Check packing ratio for objects.\n";
            exit(1);
          }
        } while (isoverlap);
        members_.push_back(member);
      }
    }
  } else if (insertion_type.compare("fill") == 0) {
    double rlength    = node["br_rod"]["rod"]["length"].as<double>();
    double max_length = node["br_rod"]["rod"]["max_length"].as<double>();
    double min_length = node["br_rod"]["rod"]["min_length"].as<double>();
    double diameter   = node["br_rod"]["rod"]["diameter"].as<double>();

    std::cout << std::setw(25) << std::left << "   length:" << std::setw(10)
      << std::left << rlength << std::endl;
    std::cout << std::setw(25) << std::left << "   max length:" << std::setw(10)
      << std::left << max_length << std::endl;
    std::cout << std::setw(25) << std::left << "   min length:" << std::setw(10)
      << std::left << min_length << std::endl;
    std::cout << std::setw(25) << std::left << "   diameter:" << std::setw(10)
      << std::left << diameter << std::endl;

    params_->rod_length = rlength;
    params_->max_rod_length = max_length;
    params_->min_rod_length = min_length;
    max_length_ = max_length;
    min_length_ = min_length;
    params_->rod_diameter = diameter;

    // Try just inserting as many rods as we can
    int nrods = 0;
    bool inserting = true;
    while(inserting) {
      BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      member->SetColor(color, draw_type);
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
            /*if (debug_trace) {
              printf("Overlap detected [oid: %d,%d], [cid: %d, %d] -> (%2.2f < %2.2f)\n",
                      part1->GetOID(), part2->GetOID(), part1->GetCID(), part2->GetCID(), idm.dr_mag2, diameter2);
            }*/
            // We can just call init again to get new random numbers
            member->Init();
          }
        } // check against current members
        if (numoverlaps > params_->max_overlap) {
          std::cout << "Done inserting!\n";
          inserting = false;
          delete member;
          break;
        }
      } while(isoverlap);
      if (!isoverlap) {
        members_.push_back(member);
        nrods++;
      }
    }
    params_->n_rod = nrods;
    std::cout << "Inserted: " << nrods << " members (" << members_.size() << ")\n";
  } else if (insertion_type.compare("oriented") == 0) {

    printf("Hey there!\n");
    //int orientation = node["br_rod"]["rod"]["orientation"].as<int>();

    int nrods         = node["br_rod"]["rod"]["num"].as<int>();
    double rlength    = node["br_rod"]["rod"]["length"].as<double>();
    double max_length = node["br_rod"]["rod"]["max_length"].as<double>();
    double diameter   = node["br_rod"]["rod"]["diameter"].as<double>();

    std::cout << std::setw(25) << std::left << "   n rods:" << std::setw(10)
      << std::left << nrods << std::endl;
    std::cout << std::setw(25) << std::left << "   length:" << std::setw(10)
      << std::left << rlength << std::endl;
    std::cout << std::setw(25) << std::left << "   max length:" << std::setw(10)
      << std::left << max_length << std::endl;
    std::cout << std::setw(25) << std::left << "   diameter:" << std::setw(10)
      << std::left << diameter << std::endl;

    params_->n_rod = nrods;
    params_->rod_length = rlength;
    params_->max_rod_length = max_length;
    max_length_ = max_length;
    params_->rod_diameter = diameter;


    for (int i = 0; i < nrods; ++i) {
      BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
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



              /*if (debug_trace) {
                printf("Overlap detected [oid: %d,%d], [cid: %d, %d] -> (%2.2f < %2.2f)\n",
                        part1->GetOID(), part2->GetOID(), part1->GetCID(), part2->GetCID(), idm.dr_mag2, diameter2);
              }*/
              // We can just call init again to get new random numbers
              member->Init();
            }
          } // check against current members
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << nrods;
            std::cout << ".  Check packing ratio for objects.\n";
            exit(1);
          }
        } while (isoverlap);
        members_.push_back(member);
      }
    }

  } else {
    printf("nope, not yet\n");
    exit(1);
  }
}

void BrRodSpecies::ConfiguratorSpindle(int ispb, int spb_oid,
                                       const double* const r_spb,
                                       const double* const u_spb,
                                       const double* const v_spb,
                                       const double* const w_spb,
                                       al_set *anchors) {
  char *filename = params_->config_file;
  //std::cout << "BrRod species\n";

  YAML::Node node = YAML::LoadFile(filename);

  // XXX FIXME default to orientation draw
  int draw_type = 1;
  double color[3] = {1.0, 0.0, 0.0};
  int nrods         = (int)node["spb"][ispb]["mt"].size();
  //std::cout << "nrods: " << nrods << std::endl;
  double max_length = node["spb"][ispb]["properties"]["mt_max_length"].as<double>();
  //std::cout << std::setprecision(16) << "max_length: " << max_length << std::endl;
  double min_length = node["spb"][ispb]["properties"]["mt_min_length"].as<double>();
  double diameter   = node["spb"][ispb]["properties"]["mt_diameter"].as<double>();

  params_->max_rod_length = max_length;
  params_->min_rod_length = min_length;
  max_length_ = max_length;
  min_length_ = min_length;
  params_->rod_diameter = diameter;

  double spb_dia = node["spb"][ispb]["properties"]["attach_diameter"].as<double>();
  double r0 = node["spb"][ispb]["properties"]["r0"].as<double>();

  double conf_diameter = space_->unit_cell[0][0];
  double conf_radius = 0.5 * conf_diameter;
  int ndim = space_->n_dim;

  double m_u_spb[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < ndim; ++i) m_u_spb[i] = r_spb[i]/conf_radius;

  std::cout << "u pointing to spb: (" << std::setprecision(16)
    << m_u_spb[0] << ", " << m_u_spb[1] << ", " << m_u_spb[2] << ")\n";

  for (int irod = 0; irod < nrods; ++irod) {

    // Break into 2 sections, one for random insertion, the other for rtp
    std::string insert_type = node["spb"][ispb]["mt"][irod]["insertion_type"].as<std::string>();
    double mlength = node["spb"][ispb]["mt"][irod]["length"].as<double>();

    double u_bond[3] = {0.0, 0.0, 0.0};
    double v_bond[3] = {0.0, 0.0, 0.0};
    double x_bond[3] = {0.0, 0.0, 0.0};

    // Create the rod
    BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());

    if (insert_type.compare("random") == 0) {
      int overlap;
      do {
        // generate random anchor position within cell
        double alpha;
        double alpha_max = asin(spb_dia/conf_diameter);
        double pos_vec[3];

        do {
          generate_random_unit_vector(ndim, pos_vec, rng_.r);
          alpha = acos(dot_product(ndim, m_u_spb, pos_vec));
        } while (alpha > alpha_max);

        // Insert with random orientation
        double rmag2;
        do {
          rmag2 = 0.0;
          generate_random_unit_vector(ndim, u_bond, rng_.r);
          // Get the - end of the bond, supposedly
          for (int i = 0; i < ndim; ++i) {
            v_bond[i] = mlength * u_bond[i];
            rmag2 += SQR(conf_radius * pos_vec[i] + v_bond[i] + u_bond[i] * r0);
          }

        } while (rmag2 >= SQR(conf_radius - 0.5) || dot_product(ndim, m_u_spb, u_bond) > 0);

        //std::cout << "found an rmag2 we liked\n";
        // Use helper functions to generate rod positions before we actually create
        // and insert the rod (like bob)
        for (int i = 0; i < ndim; ++i) {
          x_bond[i] = conf_radius * pos_vec[i] +
            0.5 * v_bond[i] +
            u_bond[i] * r0;
        }

        // Update the Rod
        member->InitConfigurator(x_bond, u_bond, mlength);
        member->SetColor(color, draw_type);

        // Check for overlaps
        overlap = 0;
        double s_bond[3] = {0.0, 0.0, 0.0};

        for (auto rodit = members_.begin(); rodit != members_.end(); ++rodit) {
          auto part1 = member->GetSimples()[0];
          auto part2 = (*rodit)->GetSimples()[0];
          interactionmindist idm;
          MinimumDistance(part1, part2, idm, space_->n_dim, space_->n_periodic, space_);
          double diameter2 = diameter*diameter;
          if (idm.dr_mag2 < diameter2) {
            overlap++;
          }
        }

        //for (int jrod = 0; jrod < irod-1; ++jrod) {
        //  auto mrod = members_[jrod];
        //  auto part1 = member->GetSimples()[0];
        //  double mrod_x[3] = {0.0, 0.0, 0.0};
        //  double mrod_s[3] = {0.0, 0.0, 0.0};
        //  double mrod_u[3] = {0.0, 0.0, 0.0};
        //  double mrod_l = mrod->GetLength();
        //  std::copy(mrod->GetPosition(), mrod->GetPosition()+3, mrod_x);
        //  std::copy(mrod->GetScaledPosition(), mrod->GetScaledPosition()+3, mrod_s);
        //  std::copy(mrod->GetOrientation(), mrod->GetOrientation()+3, mrod_u);
        //  overlap += sphero_overlap(ndim, 0, space_->unit_cell, 0.0, 1.0,
        //      x_bond,
        //      s_bond,
        //      u_bond,
        //      mlength,
        //      mrod_x,
        //      mrod_s,
        //      mrod_u,
        //      mrod_l
        //      );
        //}
      } while (overlap != 0);
    } else if (insert_type.compare("rtp") == 0) {
      // r theta phi insert
      double rod_theta = node["spb"][ispb]["mt"][irod]["theta"].as<double>();
      double rod_phi   = node["spb"][ispb]["mt"][irod]["phi"].as<double>();

      double orient[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < ndim; ++i) {
        orient[i] = node["spb"][ispb]["mt"][irod]["u"][i].as<double>();
      }
      double umag = sqrt(dot_product(ndim, orient, orient));
      for (int i = 0; i < ndim; ++i) {
        u_bond[i] = orient[i] /= umag;
      }

      x_bond[0] = conf_radius * sin(rod_theta) * cos(rod_phi) + (0.5 * mlength + r0) * u_bond[0];
      x_bond[1] = conf_radius * sin(rod_theta) * sin(rod_phi) + (0.5 * mlength + r0) * u_bond[1];
      x_bond[2] = conf_radius * cos(rod_theta) + (0.5 * mlength + r0) * u_bond[2];
    }

    // Insert the rod
    //BrRod *member = new BrRod(params_, space_, gsl_rng_get(rng_.r), GetSID());
    //member->InitConfigurator(x_bond, u_bond, mlength);
    //member->SetColor(color, draw_type);
    //member->Dump();
    member->SetAnchors(anchors);
    members_.push_back(member);

    // Add to the anchor list!!!!
    anchor_t new_anchor;
    new_anchor.idx_base_ = spb_oid;
    new_anchor.idx_other_ = member->GetSimples()[0]->GetOID();
    for (int i = 0; i < ndim; ++i) {
      new_anchor.pos0_[i] = x_bond[i] - 
        u_bond[i] * (0.5 * mlength + r0);
      new_anchor.pos1_[i] = x_bond[i] -
        mlength * u_bond[i] * 0.5;
    }
    new_anchor.pos_rel1_[0] = -0.5 * mlength;
    new_anchor.pos_rel1_[1] = new_anchor.pos_rel1_[2] = 0.0;
    double r_rel[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < ndim; ++i) {
      r_rel[i] = new_anchor.pos0_[i] - r_spb[i];
    }
    double u_proj = dot_product(ndim, r_rel, u_spb);
    double v_proj = dot_product(ndim, r_rel, v_spb);
    double w_proj = dot_product(ndim, r_rel, w_spb);
    for (int i = 0; i < 3; ++i) {
      new_anchor.pos_rel0_[i] =
        u_spb[i] * u_proj +
        v_spb[i] * v_proj +
        w_spb[i] * w_proj;
    }
    (*anchors)[spb_oid].push_back(new_anchor);

    // Tell us what just happened
    std::cout << "   New bond " << member->GetOID() << std::endl;
    std::cout << std::setprecision(16)
      << "      x (" << x_bond[0] << ", " << x_bond[1] << ", " << x_bond[2] << ")\n"
      << "      u (" << u_bond[0] << ", " << u_bond[1] << ", " << u_bond[2] << ")\n"
      << "      binding site (" << new_anchor.pos0_[0] << ", " << new_anchor.pos0_[1] << ", " << new_anchor.pos0_[2] << ")\n"
      << "      binding site rel (" << new_anchor.pos_rel0_[0] << ", " << new_anchor.pos_rel0_[1] << ", " << new_anchor.pos_rel0_[2] << ")\n"
      << "      tip site     (" << new_anchor.pos1_[0] << ", " << new_anchor.pos1_[1] << ", " << new_anchor.pos1_[2] << ")\n"
      << "      tip site rel (" << new_anchor.pos_rel1_[0] << ", " << new_anchor.pos_rel1_[1] << ", " << new_anchor.pos_rel1_[2] << ")\n"
      << "      length " << mlength << "\n"
      << "      parent anchor = " << spb_oid << std::endl;
  }
}

//Only reading out bound positions but this you might want to change
void BrRod::WritePosit(std::fstream &op){
  for (auto& velem : v_elements_)
    velem.WritePosit(op);
}

void BrRod::ReadPosit(std::fstream &ip){
  for (auto& velem : v_elements_)
    velem.ReadPosit(ip);
}

void BrRodSpecies::CreateTestRod(BrRod **rod,
                                 int ndim,
                                 std::vector<Simple*>* simples,
                                 std::unordered_map<int, int>* oid_position_map,
                                 const std::string &filename,
                                 const std::string &modulename,
                                 const std::string &unitname,
                                 const std::string &rodname,
                                 int itest) {
  YAML::Node node = YAML::LoadFile(filename);
  BrRod *mrod = *rod;
  // Load the rod
  double xr[3] = {0.0, 0.0, 0.0};
  double ur[3] = {0.0, 0.0, 0.0};
  std::ostringstream xrodname;
  xrodname << "x_" << rodname;
  std::ostringstream urodname;
  urodname << "u_" << rodname;
  std::ostringstream lrodname;
  lrodname << "l_" << rodname;
  //std::cout << "Checking location\n";
  for (int idim = 0; idim < ndim; ++idim) {
    xr[idim] = node[modulename][unitname.c_str()]["test"][itest][xrodname.str().c_str()][idim].as<double>();
    ur[idim] = node[modulename][unitname.c_str()]["test"][itest][urodname.str().c_str()][idim].as<double>();
  }
  //std::cout << "Checking length\n";
  double lrod = node[modulename][unitname.c_str()]["test"][itest][lrodname.str().c_str()].as<double>();
  //std::cout << "TEST ROD CREATE: \n";
  //std::cout << std::setprecision(16) << "x: (" << xr[0] << ", " << xr[1] << ", " << xr[2] << "), ";
  //std::cout << std::setprecision(16) << "u: (" << ur[0] << ", " << ur[1] << ", " << ur[2] << "), ";
  //std::cout << std::setprecision(16) << "l: " << lrod << std::endl;
  mrod->InitConfigurator(xr, ur, lrod);
  mrod->UpdateRodLength(0.0);
  //mrod->Dump();

  // Add the simples to the correct location and the oid position map
  std::vector<Simple*> sim_vec = mrod->GetSimples();
  for (int i = 0; i < sim_vec.size(); ++i) {
    simples->push_back(sim_vec[i]);
    (*oid_position_map)[sim_vec[i]->GetOID()] = simples->size() -1;
  }

  // Print to make sure working
  /*std::cout << "simples: \n";
  for (int i = 0; i < simples->size(); ++i) {
    std::cout << "[" << i << "] -> OID " << (*simples)[i]->GetOID();
    std::cout << " <--> " << (*oid_position_map)[(*simples)[i]->GetOID()] << std::endl;
  }*/
}

void BrRodSpecies::CreateTestRod(BrRod **rod,
                                 int ndim,
                                 std::vector<Simple*>* simples,
                                 std::unordered_map<int, int>* oid_position_map,
                                 YAML::Node *subnode) {
  YAML::Node node = *subnode;
  BrRod *mrod = *rod;
  // Load the rod
  std::cout << "Node:\n" << node << std::endl;
  double xr[3] = {0.0, 0.0, 0.0};
  double ur[3] = {0.0, 0.0, 0.0};
  for (int idim = 0; idim < ndim; ++idim) {
    xr[idim] = node["x"][idim].as<double>();
    ur[idim] = node["u"][idim].as<double>();
  }
  double lr = node["l"].as<double>();
  mrod->InitConfigurator(xr, ur, lr);
  mrod->UpdateRodLength(0.0);

  std::vector<Simple*> sim_vec = mrod->GetSimples();
  for (int i = 0; i < sim_vec.size(); ++i) {
    simples->push_back(sim_vec[i]);
    (*oid_position_map)[sim_vec[i]->GetOID()] = simples->size() -1;
  }
}
