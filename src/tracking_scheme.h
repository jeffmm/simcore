#ifndef _SIMCORE_TRACKING_SCHEME_H_
#define _SIMCORE_TRACKING_SCHEME_H_

#include "auxiliary.h"
#include "interaction.h"
#include "potential_base.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class TrackingScheme {
  public:

    TrackingScheme() {}
    virtual ~TrackingScheme() {

    }

    virtual void Init(space_struct *pSpace,
                      PotentialBase *pPotentialBase,
                      std::vector<interaction_t> *pInteractions,
                      YAML::Node *pNode);
    virtual void Print();

  protected:

    // Inputs that are needed
    int ndim_;
    int nperiodic_;
    int nthreads_;

    SID sid0_;
    SID sid1_;
    ptype type_;

    std::string name_ = "TrackingScheme";

    space_struct *space_;
    PotentialBase *pbase_;

    std::vector<interaction_t> *interactions_;
};

#endif
