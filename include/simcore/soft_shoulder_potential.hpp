#ifndef _SIMCORE_SOFT_SHOULDER_POTENTIAL_H_
#define _SIMCORE_SOFT_SHOULDER_POTENTIAL_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "potential_base.hpp"

class SoftShoulderPotential : public PotentialBase {
 protected:
  double eps_, rs_, Rs8inv_, a_;

 public:
  SoftShoulderPotential() {}
  void CalcPotential(Interaction &ix) {
    double rmag = sqrt(ix.dr_mag2);
    double r7 = pow(rmag, 7);
    double r8 = r7 * rmag;
    double R = ix.buffer_mag;
    double R8inv = 1.0 / pow(R, 8);
    double exp1 = eps_ * exp(-r8 * R8inv);
    double exp2 = eps_ * a_ * exp(-r8 * Rs8inv_);
    double ffac = -8.0 * r7 * (exp1 * R8inv + exp2 * Rs8inv_);
    double *dr = ix.dr;
    // Cut off the force at fcut
    if (ABS(ffac) > max_force_) {
      MaxForceViolation();
      ffac = SIGNOF(ffac) * max_force_;
    }
    for (int i = 0; i < n_dim_; ++i) {
      ix.force[i] = ffac * dr[i] / rmag;
    }
    for (int i = 0; i < n_dim_; ++i)
      for (int j = 0; j < n_dim_; ++j)
        ix.stress[n_dim_ * i + j] = -dr[i] * ix.force[j];
    ix.pote = exp1 + exp2;
  }

  void InitPotentialParams(system_parameters *params) {
    // Initialize potential params
    eps_ = params->ss_eps;
    a_ = params->ss_a;
    double Rs = params->ss_rs;

    // For SoftShoulderPotential potentials, the rcutoff is
    // restricted to be at 2^(1/6)sigma

    Rs8inv_ = 1.0 / pow(Rs, 8);
    rcut_ = 1.4 * Rs;  // goes quickly to zero after 1.2*rs_
    rcut2_ = rcut_ * rcut_;
  }
};

#endif
