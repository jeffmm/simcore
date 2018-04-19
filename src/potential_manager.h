#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "wca_potential.h"
#include "r2_potential.h"
#include "soft_shoulder_potential.h"
#include "soft_potential.h"

class PotentialManager {
  private:
    // Potentials
    WCAPotential wca_;
    SoftPotential soft_;
    //R2Potential r2pot_;
    //SoftShoulderPotential sspot_;
    PotentialBase * pot_;
    potential_type pot_type_;
  public:
    void InitPotentials(system_parameters *params) {
      pot_type_ = potential_type::_from_string(params->potential.c_str());
      if (pot_type_ == +potential_type::wca) {
        pot_ = &wca_;
      }
      else if (pot_type_ == +potential_type::soft) {
        pot_ = &soft_;
      }
      pot_->Init(params);
      //wca_.Init(params);
      //r2pot_.Init(params);
      //sspot_.Init(params);
    }
    void CalcPotential(Interaction *ix) {
      pot_->CalcPotential(ix);
    }
    double GetRCut2() {return pot_->GetRCut2();}
};


#endif
