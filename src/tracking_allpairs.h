#ifndef _SIMCORE_TRACKING_ALLPAIRS_H_
#define _SIMCORE_TRACKING_ALLPAIRS_H_

#include "auxiliary.h"
#include "tracking_base.h"

class TrackingAllPairs : public TrackingBase {
  public:
    TrackingAllPairs() {
      name_ = "TrackingAllPairs";
    }
    virtual ~TrackingAllPairs() {}

    virtual void CreateSubstructure(double pRcut, nl_list** pNeighbors);
    virtual void UpdateTracking(bool pForceUpdate = false);
    virtual void print();
    virtual void dump();

  protected:
    bool first_ = true;
};

#endif
