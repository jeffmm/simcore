#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "max_force_potential.hpp"
#include "r2_potential.hpp"
#include "soft_potential.hpp"
#include "soft_shoulder_potential.hpp"
#include "wca_potential.hpp"

class PotentialManager {
 private:
  // Potentials
  WCAPotential wca_;
  SoftPotential soft_;
  MaxForcePotential max_;
  // R2Potential r2pot_;
  // SoftShoulderPotential sspot_;
  PotentialBase *pot_;
  potential_type pot_type_;

 public:
  void InitPotentials(system_parameters *params) {
    /*
     * Check to see what kind of potential we are using,
     * and set that potential to be our default
     */
    pot_type_ = potential_type::_from_string(params->potential.c_str());
    if (pot_type_ == +potential_type::wca) {
      pot_ = &wca_;
      /*
       * Since WCA can result in infinite forces,
       * we initialize the max force potential in
       * case we have an overlap of objects
       */
      max_.Init(params);
    } else if (pot_type_ == +potential_type::soft) {
      pot_ = &soft_;
    }
    pot_->Init(params);
  }
  void CalcPotential(Interaction *ix) {
    /*
     *  If dr_mag2 is very very small (or zero!) we
     *  are overlapping. This can be very dangerous
     *  for infinite potentials, so enforce a cutoff.
     */
    if (ix->dr_mag2 < 1e-12 && pot_type_ == +potential_type::wca) {
      max_.CalcPotential(ix);
      return;
    }
    pot_->CalcPotential(ix);
  }
  double GetRCut2() { return pot_->GetRCut2(); }
};

#endif
