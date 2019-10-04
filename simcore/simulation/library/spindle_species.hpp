#ifndef _SIMCORE_SPINDLE_SPECIES_H_
#define _SIMCORE_SPINDLE_SPECIES_H_

#include "species.hpp"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif
#include "spindle.hpp"

class SpindleSpecies : public Species<Spindle> {
protected:
  bool midstep_;

public:
  SpindleSpecies() : Species() {
    SetSID(species_id::spindle);
    midstep_ = true;
  }
  void Init(system_parameters *params, species_parameters *sparams,
            space_struct *space) {
    Species::Init(params, sparams, space);
    sparams_ = &(params_->spindle);
  }
  void UpdatePositions() {
    for (auto it = members_.begin(); it != members_.end(); ++it) {
      it->UpdatePosition(midstep_);
    }
    midstep_ = !midstep_;
  }
};

#endif // _SIMCORE_SPINDLE_SPECIES_H_
