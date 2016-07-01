#ifndef _SIMCORE_TEST_POTENTIAL_H_
#define _SIMCORE_TEST_POTENTIAL_H_

#include "auxiliary.h"
#include "potential_base.h"

// Potential base class for all external(or maybe even internal) potentials
class TestPotential : public PotentialBase {
  public:
    TestPotential(space_struct *pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut) {
      pot_name_ = "Test Potential";
    }
    virtual ~TestPotential() {}
    void CalcPotential(double *dr,
                       double dr_mag,
                       double buffer,
                       double *fpote) {
      double rinv = 1.0/(dr_mag-buffer);
      double rinv3 = rinv*rinv*rinv;
      double fmag = -2.0*rinv3/dr_mag;
      for (int i=0; i<n_dim_; ++i) 
        fpote[i] = fmag * dr[i];
      fpote[n_dim_] = rinv*rinv;
    }

};

#endif
