#ifndef _SIMCORE_POTENTIAL_BASE_H_
#define _SIMCORE_POTENTIAL_BASE_H_

#include "auxiliary.h"
#include "interaction.h"

class PotentialBase {
  protected:
    int n_dim_;
    double fcut_,
           rcut_,
           rcut2_;
  public:
    PotentialBase() {}
    double GetRCut2() {return rcut2_;}
    virtual void CalcPotential(Interaction *ix) {}
    virtual void Init(system_parameters *params) {}
};

#endif

