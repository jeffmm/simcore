template <typename T>
double const Species<T>::GetVolume() {
  double vol = 0.0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    vol += it->GetVolume();
  }
  return vol;
}

template <typename T>
void Species<T>::Reserve() {
  members_.reserve(GetNInsert());
}

template <typename T>
void Species<T>::AddMember() {
  T newmember;
  members_.push_back(newmember);
  members_.back().SetSID(GetSID());
  members_.back().Init();
  // newmember->SetColor(sparams_->color, sparams_->draw_type);
  n_members_++;
}

template <typename T>
double const Species<T>::GetDrMax() {
  double max_dr = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    double dr = it->GetDrTot();
    if (dr > max_dr) {
      max_dr = dr;
    }
  }
  return max_dr;
}

template <typename T>
void Species<T>::ZeroDrTot() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ZeroDrTot();
  }
}

template <typename T>
void Species<T>::AddMember(T newmem) {
  members_.push_back(newmem);
  newmem.SetSID(GetSID());
  n_members_++;
}

template <typename T>
void Species<T>::PopMember() {
  members_.back().Cleanup();
  members_.pop_back();
  n_members_--;
}

template <typename T>
void Species<T>::PopAll() {
  while (n_members_ > 0) {
    PopMember();
  }
}

template <typename T>
void Species<T>::Draw(std::vector<graph_struct *> *graph_array) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Draw(graph_array);
  }
}

template <typename T>
void Species<T>::UpdatePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->UpdatePosition();
  }
}

template <typename T>
void Species<T>::GetInteractors(std::vector<Object *> *ixors) {
  //ix->clear();
  std::vector<Object *> ix;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->GetInteractors(&ix);
  }
  ixors->insert(ixors->end(), ix.begin(), ix.end());
}

template <typename T>
void Species<T>::GetLastInteractors(std::vector<Object *> *ix) {
  if (members_.size() == 0) {
    error_exit(
        "Called for last interactors of species, but species has zero "
        "members\n");
  }
  members_.back().GetInteractors(ix);
}

template <typename T>
double Species<T>::GetPotentialEnergy() {
  double pe = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    pe += it->GetPotentialEnergy();
  }
  /* The total potential energy is going to be half of the
   potential energy felt by each particle. Potential energy
   is shared, so I need to avoid double counting. */
  return 0.5 * pe;
}

template <typename T>
void Species<T>::ZeroForces() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ZeroForce();
  }
}

template <typename T>
void Species<T>::Report() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Report();
  }
}

template <typename T>
int Species<T>::GetCount() {
  int count = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    count += it->GetCount();
  }
  return count;
}

template <typename T>
void Species<T>::WritePosits() {
  int size = members_.size();
  oposit_file_.write(reinterpret_cast<char *>(&size), sizeof(size));
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->WritePosit(oposit_file_);
}

template <typename T>
void Species<T>::WriteSpecs() {
  int size = members_.size();
  ospec_file_.write(reinterpret_cast<char *>(&size), sizeof(size));
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->WriteSpec(ospec_file_);
}

template <typename T>
void Species<T>::WriteCheckpoints() {
  int size = members_.size();
  std::fstream ocheck_file(checkpoint_file_, std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    std::cout << "ERROR: Output " << checkpoint_file_ << " file did not open\n";
    exit(1);
  }
  ocheck_file.write(reinterpret_cast<char *>(&size), sizeof(size));
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->WriteCheckpoint(ocheck_file);
  ocheck_file.close();
}

template <typename T>
void Species<T>::ReadPosits() {
  if (iposit_file_.eof()) {
    printf("  EOF reached\n");
    early_exit = true;
    return;
  }
  if (!iposit_file_.is_open()) {
    printf(" ERROR: Posit file unexpectedly not open! Exiting.\n");
    early_exit = true;
    return;
  }
  int size = -1;
  T *member;
  iposit_file_.read(reinterpret_cast<char *>(&size), sizeof(size));
  // Hacky workaround FIXME
  // This prevents strange errors that occasionally crop up when reading inputs
  if (size == -1) {
    early_exit = true;
    return;
  }
  if (size != n_members_) {
    T member;
    members_.resize(size, member);
    n_members_ = size;
  }
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ReadPosit(iposit_file_);
}

template <typename T>
void Species<T>::ReadPositsFromSpecs() {
  ReadSpecs();
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->SetAvgPosition();
  }
}

template <typename T>
void Species<T>::ReadCheckpoints() {
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    std::cout << "  ERROR: Output " << checkpoint_file_
              << " file did not open\n";
    exit(1);
  }
  int size = 0;
  icheck_file.read(reinterpret_cast<char *>(&size), sizeof(size));
  T member;
  members_.resize(size, member);
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ReadCheckpoint(icheck_file);
  icheck_file.close();
}

template <typename T>
void Species<T>::ReadSpecs() {
  if (ispec_file_.eof()) {
    if (HandleEOF()) {
      printf("  Switching to new spec file\n");
      return;
    } else {
      printf("  EOF reached\n");
      early_exit = true;
      return;
    }
  }
  if (!ispec_file_.is_open()) {
    printf(" ERROR: Spec file unexpectedly not open! Exiting early.\n");
    early_exit = true;
    return;
  }
  int size = -1;
  T *member;
  ispec_file_.read(reinterpret_cast<char *>(&size), sizeof(size));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (size == -1) {
    if (HandleEOF()) {
      printf("  Switching to new spec file\n");
      return;
    } else {
      printf("  EOF reached in species\n");
      early_exit = true;
      return;
    }
  }
  if (size != n_members_) {
    T member;
    members_.resize(size, member);
    n_members_ = size;
  }
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ReadSpec(ispec_file_);
}

template <typename T>
void Species<T>::ScalePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ScalePosition();
}

template <typename T>
void Species<T>::CleanUp() {
  members_.clear();
}

template <typename T>
void Species<T>::ArrangeMembers() {
  if (GetInsertionType().compare("custom") == 0)
    CustomInsert();
  else if (GetInsertionType().compare("simple_crystal") == 0)
    CrystalArrangement();
  else if (GetInsertionType().compare("centered_oriented") == 0)
    CenteredOrientedArrangement();
  else
    warning(
        "Arrangement not recognized and ArrangeMembers not overwritten by "
        "species!\n");
}

template <typename T>
void Species<T>::SetLastMemberPosition(double const *const pos) {
  members_.back().SetPosition(pos);
  members_.back().UpdatePeriodic();
}

template <typename T>
double Species<T>::GetSpecLength() {
  if (members_.size() == 0)
    return 0;
  else
    return members_[0].GetLength();
}

template <typename T>
double Species<T>::GetSpecDiameter() {
  if (members_.size() == 0)
    return 0;
  else
    return members_[0].GetDiameter();
}

template <typename T>
void Species<T>::CenteredOrientedArrangement() {
  // This is redundant for filaments, since they have already inserted
  // themselves properly.
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  u[params_->n_dim - 1] = 1.0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->SetPosition(pos);
    it->SetOrientation(u);
  }
}

template <typename T>
void Species<T>::CrystalArrangement() {
  double d = members_[0].GetDiameter();
  double l = members_[0].GetLength();
  int n_dim = params_->n_dim;
  double R = space_->radius;
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  u[n_dim - 1] = 1;
  double nX_max = floor((2.0 * R / d) - 1);
  double nY_max = floor((2.0 * R - d) / (l + d));
  double nZ_max = (n_dim == 3 ? nX_max : 1);
  double N = nX_max * nY_max * nZ_max;
  if (n_members_ > N) {
    error_exit(
        "Number of members in species exceeds maximum possible for crystal in "
        "system radius! Max possible: %d",
        (int)N);
  }
  double fraction = N / n_members_;
  if (fraction < 1) fraction = 1;
  if (l == 0) {
    if (n_dim == 2)
      nY_max = floor(pow(n_members_, 1.0 / n_dim));
    else if (n_dim == 3)
      nY_max = ceil(pow(n_members_, 1.0 / n_dim));
  }
  if (n_dim == 2) {
    nX_max = ceil(n_members_ / nY_max);
  } else if (n_dim == 3) {
    nX_max = floor(sqrt(n_members_ / nY_max));
    if (nX_max == 0) nX_max = 1;
    nZ_max = ceil(n_members_ / (nX_max * nY_max));
  }
  for (int i = 0; i < n_dim; ++i) {
    pos[i] = -R + 0.5 * d;
  }
  pos[n_dim - 1] += 0.5 * l;
  int inserted = 0;
  int insert_x = 0;
  int insert_y = 0;
  int insert_z = 0;
  bool shift = false;
  double diff_y = (2 * R - nY_max * (l + d)) / nY_max;
  double diff_x = (2 * R - nX_max * d) / nX_max;
  double diff_z = (2 * R - nZ_max * d) / nZ_max;
  int u0 = 1;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    // Check crystal orientation type
    if (params_->uniform_crystal == 0) {
      // Random orientations
      u[n_dim - 1] = (gsl_rng_uniform_int(rng_.r, 2) == 0 ? 1 : -1);
    } else if (params_->uniform_crystal == 2) {
      // Stagger orientations
      u[n_dim - 1] = -u[n_dim - 1];
    }
    // Insert object
    it->InsertAt(pos, u);
    inserted++;
    // Update next position
    pos[n_dim - 1] += l + d + diff_y;
    if (++insert_y == nY_max) {
      if (params_->uniform_crystal == 2 && shift) {
        // Stagger orientations row-wise
        u[n_dim - 1] = u0;
        u0 = -u0;
      }
      // Stagger positions column-wise
      shift = !shift;
      insert_y = 0;
      pos[n_dim - 1] =
          -R + 0.5 * (l + d) + (shift ? 0.5 * (l + d + diff_y) : 0);
      pos[0] += d + diff_x;
      if (++insert_x == nX_max) {
        if (inserted < n_members_ && n_dim == 2) {
          error_exit(
              "Ran out of room while arranging crystal in 2D! Arranged %d/%d",
              inserted, n_members_);
        } else if (n_dim == 2)
          continue;
        insert_x = 0;
        pos[0] = -R + 0.5 * d;
        pos[1] += d + diff_z;
        if (++insert_z == nZ_max && inserted < n_members_) {
          error_exit(
              "Ran out of room while arranging crystal in 3D! Arranged %d/%d",
              inserted, n_members_);
        }
      }
    }
  }
}

template <typename T>
void Species<T>::CustomInsert() {
  YAML::Node inode;
  try {
    inode = YAML::LoadFile(sparams_->insert_file);
  } catch (...) {
    std::cout << "Failed to load custom insert file " << sparams_->insert_file
              << " for species " << spec_name_ << "\n";
    error_exit("");
  }
  if (!inode[spec_name_] || inode[spec_name_].size() != n_members_) {
    std::cout << "Custom insert file for species " << spec_name_
              << " was invalid: \n";
    if (!inode[spec_name_]) {
      std::cout << "  Species ID header was missing from file\n";
    } else if (inode[spec_name_].size() != n_members_) {
      std::cout << "  There were " << inode[spec_name_].size()
                << " positions specified for " << n_members_
                << " species members\n";
    }
    error_exit("");
  }
  for (YAML::const_iterator it = inode.begin(); it != inode.end(); ++it) {
    if (it->first.as<std::string>().compare(spec_name_) != 0) continue;
    if (!it->second.IsSequence()) {
      error_exit("Custom insert file positions not specified as sequence");
    }
    int i_member = 0;
    for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
         ++jt) {
      double pos[3], u[3];
      if (!jt->IsSequence()) {
        error_exit("Custom insert position not specified as sequence");
      }
      if (jt->size() != 2 || (*jt)[0].size() != 3 || (*jt)[1].size() != 3) {
        error_exit(
            "Custom insert position not in format "
            "[[pos_x,pos_y,pos_z],[u_x,u_y,u_z]]");
      }
      for (int i = 0; i < 3; ++i) {
        pos[i] = (*jt)[0][i].as<double>();
        u[i] = (*jt)[1][i].as<double>();
      }
      members_[i_member++].InsertAt(pos, u);
    }
  }
}
