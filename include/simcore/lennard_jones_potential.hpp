#ifndef _SIMCORE_LENNARD_JONES_POTENTIAL_H_
#define _SIMCORE_LENNARD_JONES_POTENTIAL_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "potential_base.hpp"

class LennardJonesPotential : public PotentialBase {
protected:
  double eps_, sigma_, c12_, c6_;

public:
  LennardJonesPotential() {}
  void CalcPotential(Interaction &ix) {
    double rmag = sqrt(ix.dr_mag2);
    double *dr = ix.dr;
    double rinv = 1.0 / (rmag);
    double rinv2 = rinv * rinv;
    double r6 = rinv2 * rinv2 * rinv2;
    double ffac = -(12.0 * c12_ * r6 - 6.0 * c6_) * r6 * rinv;
    if (ABS(ffac) > max_force_) {
      MaxForceViolation();
    }
    for (int i = 0; i < n_dim_; ++i) {
      ix.force[i] = ffac * dr[i] * rinv;
    }
    for (int i = 0; i < n_dim_; ++i)
      for (int j = 0; j < n_dim_; ++j)
        ix.stress[n_dim_ * i + j] = -dr[i] * ix.force[j];
    ix.pote = r6 * (c12_ * r6 - c6_) + eps_;
  }
  void InitPotentialParams(system_parameters *params) {
    // Initialize potential params
    eps_ = params->wca_eps;
    sigma_ = params->wca_sig;

    // For WCAPotential potentials, the rcutoff is
    // restricted to be at 2^(1/6)sigma

    rcut_ = 2.5 * sigma_;
    rcut2_ = rcut_ * rcut_;
    c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
    c6_ = 4.0 * eps_ * pow(sigma_, 6.0);
  }
};

#endif
