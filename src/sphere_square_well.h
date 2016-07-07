#ifndef _SIMCORE_SPHERE_SQUARE_WELL_H_
#define _SIMCORE_SPHERE_SQUARE_WELL_H_

#include "auxiliary.h"
#include "potential_base.h"

class SphereSquareWell : public PotentialBase {
  protected:
    double depth_;
  public:
    SphereSquareWell() : PotentialBase(nullptr, 0.0, 0.0), depth_(0.0) {
      pot_name_ = "SphereSquareWell";
    }
    SphereSquareWell(double pDepth, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut), depth_(pDepth) {
      pot_name_ = "SphereSquareWell";
    }
    virtual void Print() {
      PotentialBase::Print();
      std::cout << "\t{depth:" << depth_ << "}\n";
    }
    // X and Y are real positions, XS and YS are scaled positions
    virtual void CalcPotential(double *dr,
                               double dr_mag,
                               double buffer,
                               double *fpote) {
      std::fill(fpote, fpote + n_dim_ + 1, 0.0);
      // 0 force, return depth when inside well
      if (dr_mag < rcut_) {
        fpote[n_dim_] = depth_;
      } else {
        fpote[n_dim_] = 0.0;
      }
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
      PotentialBase::Init(pSpace, ipot, node);

      // Now, let's look at the particular yaml node we are supposed to be interested in
      depth_  = node["potentials"][ipot]["depth"].as<double>();
      rcut_   = node["potentials"][ipot]["rcut"].as<double>();
      fcut_   = node["potentials"][ipot]["fcut"].as<double>();

      rcut2_ = rcut_*rcut_;
    }
};

#endif

