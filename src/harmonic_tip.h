#ifndef _SIMCORE_XLINK_HARMONIC_TIP_H_
#define _SIMCORE_XLINK_HARMONIC_TIP_H_

#include "auxiliary.h"
#include "potential_base.h"

class HarmonicTip : public PotentialBase {
  protected:
    double k_;
    double r_equil_;
    int tip_;
  public:
    HarmonicTip() : PotentialBase(nullptr, 0.0, 0.0), k_(0.0), r_equil_(0.0) {
      pot_name_ = "HarmonicTip";
    }
    HarmonicTip(double pK, double pRequil, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut), k_(pK), r_equil_(pRequil) {
      pot_name_ = "HarmonicTip";
    }
    virtual void Print();
    virtual void CalcPotential(interactionmindist *idm,
                               Simple *part1,
                               Simple *part2,
                               double *fpote);
    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node);

    // Specific functions
    const double GetK() {return k_;}
    const double GetRequil() {return r_equil_;}
    void SetK(const double k) {k_=k;}
    void SetRequil(const double requil) {r_equil_=requil;}
};

#endif

