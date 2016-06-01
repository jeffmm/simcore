#ifndef _CYTOSCORE_LJ126_H_
#define _CYTOSCORE_LJ126_H_

#include "auxiliary.h"
#include "potential_base.h"

class LJ126 : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
  public:
    LJ126(double pEps, double pSigma, space_struct* pSpace, double pRcut) : PotentialBase(pSpace, pRcut), eps_(pEps), sigma_(pSigma) {
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
    virtual void Print() {
        std::cout << "Lennard-Jones 12-6 potential:\n";
        PotentialBase::Print();
        std::cout << "\t{eps:" << eps_ << "}, {sigma:" << sigma_ << "}, {c6:" << c6_ << "}, {c12:" << c12_ << "}\n";
    }
    // X and Y are real positions, XS and YS are scaled positions
    virtual void CalcPotential(double *dr,
                               double dr_mag,
                               double buffer,
                               double *fpote) {
      std::fill(fpote, fpote + n_dim_ + 1, 0.0);
      double rmag = dr_mag;
      double ffac, r6, rinv;

      rinv = 1.0/(rmag*rmag);
      r6 = rinv*rinv*rinv;

      ffac = -(12.0*c12_*r6 - 6.0*c6_)*r6/rmag/dr_mag;
      for (int i = 0; i < n_dim_; ++i) 
        fpote[i] = dr[i]*ffac;
      fpote[n_dim_] = r6*(c12_*r6 - c6_);
    }
};

#endif

