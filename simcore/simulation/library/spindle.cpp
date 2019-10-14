#include "spindle.hpp"

Spindle::Spindle(unsigned long seed) : BrRod(seed) {
  SetSID(species_id::spindle);
}

void Spindle::SetParameters() {
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  length_ = sparams_->length;
  n_filaments_bud_ = sparams_->n_filaments_bud;
  n_filaments_mother_ = sparams_->n_filaments_mother;
  n_filaments_ = n_filaments_bud_ + n_filaments_mother_;
  k_spring_ = sparams_->k_spring;
  k_align_ = sparams_->k_align;
  spring_length_ = sparams_->spring_length;
  alignment_potential_ = (sparams_->alignment_potential ? true : false);
  spb_diameter_ = sparams_->spb_diameter;
  // Force spherocylinder insertion options to match spindle insertion options
  // FIXME
  SetDiffusion();
}

void Spindle::InitFilamentParameters(filament_parameters *fparams) {
  fparams_ = fparams;
  anchor_distance_ = 0.5 * diameter_ + (0.5 + 1e-3) * fparams_->diameter;
  InsertRod(sparams_->insertion_type,
            length_ + anchor_distance_ + 2 * fparams_->min_bond_length);
  UpdatePeriodic();
  GetBodyFrame();
  GenerateNucleationSites();

  for (int i_fil = 0; i_fil < n_filaments_; ++i_fil) {
    InsertFilament(i_fil);
  }
}

void Spindle::Init(spindle_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  filaments_.reserve(n_filaments_);
  nuc_sites_.reserve(n_filaments_);
  fil_sites_.reserve(n_filaments_);
}

// Returns true if successful, false otherwise
void Spindle::InsertFilament(int i_fil) {
  Filament fil(rng_.GetSeed());
  filaments_.push_back(fil);
  filaments_.back().Init(fparams_);
  int sign = (i_fil < n_filaments_bud_ ? 1 : -1);
  const double *const u = nuc_sites_[i_fil].GetOrientation();
  const double *const site_pos = nuc_sites_[i_fil].GetPosition();
  double new_pos[3] = {0, 0, 0};
  for (int i = 0; i < n_dim_; ++i) {
    new_pos[i] = site_pos[i] + fparams_->min_bond_length * u[i];
  }
  filaments_.back().InsertAt(new_pos, u);
  fil_sites_.push_back(filaments_.back().GetSite(0));
}

int Spindle::GetCount() {
  int count = 1;
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    count += it->GetCount();
  }
  return count;
}

void Spindle::ZeroForce() {
  Object::ZeroForce();
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    it->ZeroForce();
  }
  for (auto it = nuc_sites_.begin(); it != nuc_sites_.end(); ++it) {
    it->ZeroForce();
  }
}

void Spindle::GenerateNucleationSites() {
  for (int i = 0; i < n_filaments_; ++i) {
    Site s(rng_.GetSeed());
    nuc_sites_.push_back(s);
    nuc_sites_.back().SetDiameter(1); // For graphing

    double theta =
        2.0 * (rng_.RandomUniform() - 0.5) * atan(spb_diameter_ / diameter_);
    double phi =
        2.0 * (rng_.RandomUniform() - 0.5) * atan(spb_diameter_ / diameter_);
    nuc_sites_.back().SetTheta(theta);
    nuc_sites_.back().SetPhi(phi);
    Logger::Trace("Theta: %2.2f, phi: %2.2f", theta, phi);
  }
  ResetSitePositions();
}

void Spindle::ResetSitePositions() {
  double neg_orientation[3];
  for (int i = 0; i < 3; ++i) {
    neg_orientation[i] = -orientation_[i];
  }
  for (int i_fil = 0; i_fil < n_filaments_; ++i_fil) {
    double u[3];
    if (i_fil >= n_filaments_bud_) {
      std::copy(neg_orientation, neg_orientation + 3, u);
    } else {
      std::copy(orientation_, orientation_ + 3, u);
    }
    if (n_dim_ == 2) {
      rotate_vector(u, &body_frame_[0], nuc_sites_[i_fil].GetTheta(), n_dim_);
    } else {
      rotate_vector(u, &body_frame_[0], nuc_sites_[i_fil].GetTheta(), n_dim_);
      rotate_vector(u, &body_frame_[3], nuc_sites_[i_fil].GetPhi(), n_dim_);
    }
    normalize_vector(u, n_dim_);
    int sign = (i_fil >= n_filaments_bud_ ? -1 : 1);
    double new_pos[3] = {0, 0, 0};
    for (int i = 0; i < n_dim_; ++i) {
      new_pos[i] = position_[i] + sign * 0.5 * orientation_[i] * length_ +
                   anchor_distance_ * u[i];
    }
    Logger::Trace("Initializing nucleation site at position [%2.2f %2.2f %2.2f]"
                  "and orientation [%2.2f %2.2f %2.2f]",
                  new_pos[0], new_pos[1], new_pos[2], u[0], u[1], u[2]);
    nuc_sites_[i_fil].SetOrientation(u);
    nuc_sites_[i_fil].SetPosition(new_pos);
    // nuc_sites_[i_fil].UpdatePeriodic();
  }
}

void Spindle::UpdatePosition(bool midstep) {
  midstep_ = midstep;
  ApplyForcesTorques();
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Filament>::iterator,
                        std::vector<Filament>::iterator>>
      chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = filaments_.size() / max_threads;
  filament_iterator cur_iter = filaments_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, filaments_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
        it->UpdatePosition(midstep);
      }
    }
  }
#else
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    it->UpdatePosition(midstep);
  }
#endif
  SetPrevPosition(position_);
  if (!midstep_) {
    Integrate();
  }
}

void Spindle::ApplyForcesTorques() {
  ApplyNucleationSiteForces();
  if (!midstep_) {
    ApplySpindleForces();
  }
}

void Spindle::ApplyNucleationSiteForces() {
  for (int i_fil = 0; i_fil < n_filaments_; ++i_fil) {
    const double *const force = nuc_sites_[i_fil].GetForce();
    const double *const site_pos = nuc_sites_[i_fil].GetPosition();
    const double *const fil_pos = fil_sites_[i_fil]->GetPosition();
    double fil_forces[3] = {0, 0, 0};
    for (int i = 0; i < n_dim_; ++i) {
      fil_forces[i] = -k_spring_ * (site_pos[i] - fil_pos[i]);
    }
    fil_sites_[i_fil]->SubForce(fil_forces);
    nuc_sites_[i_fil].AddForce(fil_forces);
  }
}

void Spindle::ApplySpindleForces() {
  for (auto it = nuc_sites_.begin(); it != nuc_sites_.end(); ++it) {
    // First calculate translational forces
    double dr[3] = {0, 0, 0};
    const double *const site_pos = it->GetPosition();
    const double *const site_force = it->GetForce();
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] += site_force[i];
      dr[i] = site_pos[i] - position_[i];
    }
    // Then calculate torques (always use 3 dimensions for torques)
    double i_torque[3] = {0, 0, 0};
    cross_product(dr, site_force, i_torque, 3);
    for (int i = 0; i < 3; ++i) {
      torque_[i] += i_torque[i];
    }
  }
}

void Spindle::Integrate() {
  BrRod::Integrate();
  // Now I need to update the positions of the anchors based on the new
  // positions and orientations using their euler angles
  ResetSitePositions();
}

void Spindle::GetInteractors(std::vector<Object *> &ix) {
  ix.push_back(this);
  for (auto it = filaments_.begin(); it != filaments_.end(); ++it) {
    it->GetInteractors(ix);
  }
}

void Spindle::Draw(std::vector<graph_struct *> &graph_array) {
  Object::Draw(graph_array);
  for (auto fil = filaments_.begin(); fil != filaments_.end(); ++fil) {
    fil->Draw(graph_array);
  }
  for (auto it = nuc_sites_.begin(); it != nuc_sites_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void Spindle::UpdateDrTot() {
  if (midstep_) {
    return;
  }
  Object::UpdateDrTot();
  for (auto it = filaments_.begin(); it != filaments_.end(); ++it) {
    const double dr = it->GetDrTot();
    if (dr > dr_tot_) {
      dr_tot_ = dr;
    }
  }
}

void Spindle::ZeroDrTot() {
  Object::ZeroDrTot();
  for (auto it = filaments_.begin(); it != filaments_.end(); ++it) {
    it->ZeroDrTot();
  }
}

const double Spindle::GetDrTot() {
  UpdateDrTot();
  return dr_tot_;
}

const bool Spindle::CheckInteractorUpdate() {
  for (auto it = filaments_.begin(); it != filaments_.end(); ++it) {
    if (it->CheckInteractorUpdate()) {
      interactor_update_ = true;
    }
  }
  return Object::CheckInteractorUpdate();
}
