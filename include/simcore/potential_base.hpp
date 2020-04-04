#ifndef _SIMCORE_POTENTIAL_BASE_H_
#define _SIMCORE_POTENTIAL_BASE_H_

#include "auxiliary.hpp"
#include "interaction.hpp"

class PotentialBase {
protected:
  int n_dim_;
  double fcut_, rcut_, rcut2_;
  bool fcut_violation_ = false;

public:
  PotentialBase() {}
  double GetRCut2() { return rcut2_; }
  bool CheckFcutViolation() { return fcut_violation_; }
  virtual void CalcPotential(Interaction &ix) {}
  virtual void Init(system_parameters *params) {}
  virtual void ResetFcutViolation() { fcut_violation_ = false; }
};

#endif
