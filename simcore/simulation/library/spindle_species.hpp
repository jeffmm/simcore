#ifndef _SIMCORE_SPINDLE_SPECIES_H_
#define _SIMCORE_SPINDLE_SPECIES_H_

#include "species.hpp"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif
#include "spindle.hpp"

class SpindleSpecies : public Species<Spindle, species_id::spindle> {
protected:
  bool midstep_;

public:
  SpindleSpecies() : Species() {
    midstep_ = true;
  }
  void Init(system_parameters *params, species_base_parameters *sparams,
            space_struct *space);
  void UpdatePositions();
};

#endif // _SIMCORE_SPINDLE_SPECIES_H_
