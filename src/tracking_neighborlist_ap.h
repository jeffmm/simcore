#ifndef _SIMCORE_TRACKING_NEIGHBORLIST_AP_H_
#define _SIMCORE_TRACKING_NEIGHBORLIST_AP_H_

#include "auxiliary.h"
#include "tracking_base.h"

class TrackingNeighborListAP : public TrackingBase {
  public:
    TrackingNeighborListAP() {
      name_ = "TrackingNeighborListAP";
    }
    virtual ~TrackingNeighborListAP() {}

    virtual void CreateSubstructure(double pRcut, nl_list** pNeighbors);
    virtual void UpdateTracking(bool pForceUpdate = false);
    virtual void print();
    virtual void dump();

    void UpdateNeighborList();
    void AllPairsUpdate();

  protected:
    
    // Neighbor list needs...
    bool nl_update_;
    double rcut2_;
    double skin2_;
    double rcs2_;
    double half_skin2_;
};

#endif
