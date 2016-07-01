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
    WCA() : PotentialBase(nullptr, 0.0, 0.0), eps_(0.0), sigma_(0.0) {
      pot_name_ = "WCA";
    }
    WCA(double pEps, double pSigma, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut), eps_(pEps), sigma_(pSigma) {
      pot_name_ = "WCA";
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
      // Cut off the force at fcut
      if (ABS(ffac) > fcut_) {
        ffac = SIGNOF(ffac) * fcut_;
      }
      for (int i = 0; i < n_dim_; ++i) 
        fpote[i] = ffac*dr[i]/dr_mag;
      fpote[n_dim_] = r6*(c12_*r6 - c6_) + eps_;
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
        PotentialBase::Init(pSpace, ipot, node);

        // Now, let's look at the particular yaml node we are supposed to be interested in
        eps_    = node["potentials"][ipot]["eps"].as<double>();
        sigma_  = node["potentials"][ipot]["sigma"].as<double>();
        fcut_ = node["potentials"][ipot]["fcut"].as<double>();

        // For WCA potentials, the rcutoff is actually important, as it must be
        // restricted to be at 2^(1/6)sigma

        rcut_ = pow(2.0, 1.0/6.0)*sigma_;

        rcut2_ = rcut_*rcut_;
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }

};

#endif

