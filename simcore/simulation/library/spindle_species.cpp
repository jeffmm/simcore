#include "spindle_species.hpp"

void SpindleSpecies::Init(system_parameters *params, species_base_parameters *sparams,
          space_struct *space) {
  Species::Init(params, sparams, space);
}

void SpindleSpecies::UpdatePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->UpdatePosition(midstep_);
  }
  midstep_ = !midstep_;
}

