#ifndef _SIMCORE_LJ126_H_
#define _SIMCORE_LJ126_H_

#include "auxiliary.h"
#include "potential_base.h"

class LJ126 : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
    double shift_;
  public:
    LJ126() : PotentialBase(nullptr, 0.0), eps_(0.0), sigma_(0.0) {
      pot_name_ = "Lennard Jones 12-6";
    }
    LJ126(double pEps, double pSigma, space_struct* pSpace, double pRcut) : PotentialBase(pSpace, pRcut), eps_(pEps), sigma_(pSigma) {
      pot_name_ = "Lennard Jones 12-6";
      c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
      c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
      // Shift potential so it goes to zero at rcut
      double rcutinv = 1.0/pRcut;
      double rcutinv6 = pow(rcutinv, 6.0);
      shift_ = rcutinv6*(c12_*rcutinv6 - c6_);
    }
    virtual void Print() {
        PotentialBase::Print();
        std::cout << "\t{eps:" << eps_ << "}, {sigma:" << sigma_ << "}, {c6:" << c6_ << "}, {c12:" 
                  << c12_ << "}, {shift:" << shift_ << "}\n";
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
        fpote[i] = dr[i]*ffac;
      fpote[n_dim_] = r6*(c12_*r6 - c6_) - shift_;
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
        PotentialBase::Init(pSpace, ipot, node);

        // Now, let's look at the particular yaml node we are supposed to be interested in
        rcut_   = node["potentials"][ipot]["rcut"].as<double>();
        eps_    = node["potentials"][ipot]["eps"].as<double>();
        sigma_  = node["potentials"][ipot]["sigma"].as<double>();

        rcut2_ = rcut_*rcut_;
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
        // Shift potential so it goes to zero at rcut
        double rcutinv = 1.0/rcut_;
        double rcutinv6 = pow(rcutinv, 6.0);
        shift_ = rcutinv6*(c12_*rcutinv6 - c6_);
    }
};

#endif

