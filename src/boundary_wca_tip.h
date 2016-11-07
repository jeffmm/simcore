#ifndef _SIMCORE_BOUNDARY_TIP_WCA_H_
#define _SIMCORE_BOUNDARY_TIP_WCA_H_

#include "auxiliary.h"
#include "potential_base.h"

class BoundaryWCATip : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
    double shift_;
    double conf_radius_;
    double conf_radius2_;
    int tip_;
  public:
    BoundaryWCATip() : PotentialBase(nullptr, 0.0, 0.0), eps_(0.0), sigma_(0.0), conf_radius_(0.0) {
      pot_name_ = "BoundaryWCATip";
    }
    BoundaryWCATip(double pEps, double pSigma, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut), eps_(pEps), sigma_(pSigma) {
      pot_name_ = "BoundaryWCATip";
      c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
      c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
    virtual void Print();
    virtual void CalcPotential(interactionmindist *idm,
                               Simple *part1,
                               Simple *part2,
                               double *fpote);
    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node);
    virtual void Init(space_struct *pSpace, YAML::Node *subnode);
};

#endif

