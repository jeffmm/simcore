#ifndef _SIMCORE_POTENTIAL_BASE_H_
#define _SIMCORE_POTENTIAL_BASE_H_

#include "auxiliary.hpp"
#include "interaction.hpp"

class PotentialBase {
protected:
  static int n_dim_;
  static double max_force_;
  static bool fmax_violation_;

  double rcut_;
  double rcut2_;

public:
  PotentialBase() {}
  // Static functions
  static void SetNDim(int ndim);
  static bool CheckMaxForceViolation();
  static void MaxForceViolation();
  static void ResetMaxForceViolation();
  static void SetMaxForce(double fmax);
  static double GetMaxForce();

  void Init(system_parameters *params);
  double GetRCut2() { return rcut2_; }

  // Virtual functions
  virtual void CalcPotential(Interaction &ix) {}
  virtual void InitPotentialParams(system_parameters *params) = 0;
};
#endif
