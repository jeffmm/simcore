#ifndef _SIMCORE_MAX_FORCE_POTENTIAL_H_
#define _SIMCORE_MAX_FORCE_POTENTIAL_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "potential_base.hpp"

class MaxForcePotential : public PotentialBase {
 public:
  MaxForcePotential() {}
  void CalcPotential(Interaction *ix) {
    /* Check if we can generate a non-zero vector between
       the COMs of the two objects */
    double *dr = ix->dr;
    if (ix->contact1[0] || ix->contact1[1] || ix->contact1[2]) {
      for (int i = 0; i < n_dim_; ++i) {
        dr[i] = ix->contact1[i] - ix->contact2[i];
      }
    }
    double rmag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      rmag += dr[i] * dr[i];
    }
    rmag = sqrt(rmag);
    /* Assure that we have a nonzero distance defined
       between the objects */
    if (rmag < 1e-12) {
      // Abort
      return;
    }
    double rinv = 1.0 / (rmag);
    for (int i = 0; i < n_dim_; ++i) {
      ix->force[i] = fcut_ * dr[i] * rinv;
    }
    for (int i = 0; i < n_dim_; ++i)
      for (int j = 0; j < n_dim_; ++j)
        ix->stress[n_dim_ * i + j] = -dr[i] * ix->force[j];
    // XXX This needs a corresponding potential value
    ix->pote = 0;
  }

  void Init(system_parameters *params) {
    // Initialize potential params
    n_dim_ = params->n_dim;
    fcut_ = params->f_cutoff;
  }
};

#endif
