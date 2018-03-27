#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "wca_potential.h"
#include "r2_potential.h"
#include "soft_shoulder_potential.h"

class PotentialManager {
  public:
    // Potentials
    WCAPotential wca_;
    R2Potential r2pot_;
    SoftShoulderPotential sspot_;

    void InitPotentials(system_parameters *params) {
      wca_.Init(params);
      r2pot_.Init(params);
      sspot_.Init(params);
    }
};


#endif
