#ifndef _SIMCORE_HARMONIC_H_
#define _SIMCORE_HARMONIC_H_

#include "auxiliary.h"
#include "potential_base.h"

#include <iomanip>

class Harmonic : public PotentialBase {
  protected:
    double k_;
    double r_equil_;
  public:
    Harmonic() : PotentialBase(nullptr, 0.0, 0.0), k_(0.0), r_equil_(0.0) {
      pot_name_ = "Harmonic";
    }
    Harmonic(double pK, double pRequil, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut), k_(pK), r_equil_(pRequil) {
      pot_name_ = "Harmonic";
    }
    virtual void Print() {
        PotentialBase::Print();
        std::cout << "\tk:                  " << std::setprecision(16) << k_ << std::endl;
        std::cout << "\tequilibrium_length: " << std::setprecision(16) << r_equil_ << std::endl;
    }

    virtual void CalcPotential(interactionmindist *idm,
                               Simple *part1,
                               Simple *part2,
                               double *fpote) {
      std::fill(fpote, fpote + n_dim_ + 1, 0.0);
      double rmag = idm->dr_mag;
      double *dr = idm->dr;

      double k = k_;
      double ffac;
      if (r_equil_ == 0.0) {
        ffac = k;
      } else {
        ffac = k * (1.0 - r_equil_ / sqrt(dot_product(n_dim_, dr, dr)));
      }
      double u = 0.0;

      if (ABS(ffac) > fcut_)
        ffac = SIGNOF(ffac) * fcut_;
      for (int i = 0; i < n_dim_; ++i)  {
        fpote[i] = dr[i]*ffac;
        u += 0.5 * SQR(fpote[i]);
      }
      u *= 0.5 / k;
      fpote[n_dim_] = u;
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
      std::cout << "harmoic init\n";
        PotentialBase::Init(pSpace, ipot, node);

        std::cout << "specifices\n";
        // Now, let's look at the particular yaml node we are supposed to be interested in
        rcut_     = node["potentials"][ipot]["rcut"].as<double>();
        k_        = node["potentials"][ipot]["k"].as<double>();
        r_equil_  = node["potentials"][ipot]["equilibrium_length"].as<double>();

        rcut2_ = rcut_*rcut_;
    }
};

#endif

