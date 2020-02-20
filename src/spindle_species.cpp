#include <simcore/spindle_species.hpp>

SpindleSpecies::SpindleSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::spindle);
}
void SpindleSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (sparams_.num > 0 &&
      (sparams_.n_filaments_bud > 0 || sparams_.n_filaments_mother > 0)) {
    std::string fil_spec_name("spindle");
    species_base_parameters *sparams =
        parser.GetNewSpeciesParameters(species_id::filament, fil_spec_name);
    fparams_ = *dynamic_cast<filament_parameters *>(sparams);
    delete sparams;
  }
}

void SpindleSpecies::UpdatePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->UpdatePosition(midstep_);
  }
  midstep_ = !midstep_;
}

void SpindleSpecies::AddMember() {
  Species::AddMember();
  members_.back().InitFilamentParameters(&fparams_);
}

const double SpindleSpecies::GetSpecLength() const { 
  if (fparams_.dynamic_instability_flag) {
    return 2 * fparams_.min_bond_length;
  } else {
    return 1.5 * fparams_.min_bond_length;
  }
}

void SpindleSpecies::ReadSpecs() {
  if (ispec_file_.eof()) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (!ispec_file_.is_open()) {
    Logger::Warning("ERROR. Spec file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  n_members_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_members_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_members_ == -1) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (n_members_ == 0) {
    members_.clear();
  } else if (n_members_ != members_.size()) {
    Spindle s(rng_.GetSeed());
    s.Init(&sparams_);
    s.InitFilamentParameters(&fparams_);
    s.SetSID(GetSID());
    members_.resize(n_members_, s);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadSpec(ispec_file_);
  }
}



