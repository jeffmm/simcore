#ifndef _SIMCORE_WCA_H_
#define _SIMCORE_WCA_H_

#include "auxiliary.h"
#include "interaction.h"
//#include "potential_base.h"

class WCA {
  protected:
    int n_dim_;
    double eps_,
           sigma_,
           fcut_,
           c12_,
           c6_,
           rcut_,
           rcut2_,
           shift_;
  public:
    WCA() {}
    double GetRCut2() {return rcut2_;}
    void CalcPotential(Interaction *ix) {
      double rmag = sqrt(ix->dr_mag2);
      double *dr = ix->dr;
      double rinv = 1.0/(rmag);
      double rinv2 = rinv*rinv;
      double r6 = rinv2*rinv2*rinv2;
      double ffac = -(12.0*c12_*r6 - 6.0*c6_)*r6*rinv;
      // Cut off the force at fcut
      if (ABS(ffac) > fcut_) {
        ffac = SIGNOF(ffac) * fcut_;
      }
      ix->virial = 0;
      for (int i = 0; i < n_dim_; ++i) {
        ix->force[i] = ffac*dr[i]/rmag;
        ix->virial += ix->force[i]*dr[i];
      }
      ix->pote = r6*(c12_*r6 - c6_) + eps_;
    }

    void Init(system_parameters *params) {
      // Initialize potential params
      n_dim_ = params->n_dim;
      eps_    = params->wca_eps;
      sigma_  = params->wca_sig;
      fcut_ = params->f_cutoff;

      // For WCA potentials, the rcutoff is
      // restricted to be at 2^(1/6)sigma

      rcut_ = pow(2.0, 1.0/6.0)*sigma_;
      rcut2_ = rcut_*rcut_;
      c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
      c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
};

#endif

