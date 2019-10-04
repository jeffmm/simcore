#ifndef _SIMCORE_CROSSLINK_FACTORY_H_
#define _SIMCORE_CROSSLINK_FACTORY_H_

#include "crosslink_species.hpp"

class CrosslinkFactory {
  public:
    CrosslinkSpecies* CreateCrosslinkSpecies(const species_id sid) const {
      if (sid == +species_id::crosslink) {
        return new CrosslinkSpecies();
      }
      return nullptr;
    }
};

#endif
