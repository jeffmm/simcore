#ifndef _SIMCORE_WCA_H_
#define _SIMCORE_WCA_H_

#include "auxiliary.h"
#include "potential_base.h"

class WCA : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
    double shift_;
  public:
    WCA(double pEps, double pSigma, space_struct* pSpace, double pRcut) : PotentialBase(pSpace, pRcut), eps_(pEps), sigma_(pSigma) {
      pot_name_ = "Lennard Jones 12-6";
      c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
      c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
    virtual void Print() {
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

      ffac = -(12.0*c12_*r6 - 6.0*c6_)*r6*rinv;
      for (int i = 0; i < n_dim_; ++i) 
        fpote[i] = ffac*dr[i]/dr_mag;
      fpote[n_dim_] = r6*(c12_*r6 - c6_) + eps_;
    }
};

#endif

