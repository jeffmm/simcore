#ifndef _SIMCORE_SPINDLE_SPECIES_H_
#define _SIMCORE_SPINDLE_SPECIES_H_

#include "spindle.hpp"
#include "species.hpp"

class SpindleSpecies : public Species<Spindle, species_id::spindle> {
protected:
  bool midstep_ = true;
  filament_parameters fparams_;

public:
  SpindleSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void UpdatePositions();
  void AddMember();
  const double GetSpecLength() const;
};



#endif // _SIMCORE_SPINDLE_SPECIES_H_
