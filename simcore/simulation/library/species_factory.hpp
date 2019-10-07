#ifndef _SIMCORE_SPECIES_FACTORY_H_
#define _SIMCORE_SPECIES_FACTORY_H_

#include "filament_species.hpp"
#include "crosslink_species.hpp"
//#include "br_bead_species.hpp"
//#include "bead_spring_species.hpp"
//#include "spindle_species.hpp"
//#include "spherocylinder_species.hpp"

class SpeciesFactory {
  public:
    SpeciesBase* CreateSpecies(const species_id sid) const {
      if (sid == +species_id::filament) {
        return new FilamentSpecies();
      //} else if (sid == +species_id::br_bead) {
        //return new BrBeadSpecies();
      } else if (sid == +species_id::crosslink) {
        return new CrosslinkSpecies();
      }
      Logger::Error("Species ID not recognized in SpeciesFactory!");
      return nullptr;
    }
};

#endif
