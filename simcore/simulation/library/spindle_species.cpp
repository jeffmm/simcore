#include "spindle_species.hpp"

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



