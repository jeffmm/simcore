#ifndef _SIMCORE_BOND_H_
#define _SIMCORE_BOND_H_ 

#include "site.h"

// Bonds, ie graph edges
class Bond : public Object {
  protected:
    Site * sites_[2];
    double dr_[3],
           equil_length_,
           k_spring_,
           body_frame_[6];
  public:
    Bond();
    void Init(Site * s1, Site * s2);
    void ReInit();
    void Report();
    void ReportSites();
    //void GetBodyFrame();
    //void UpdateOrientation();
    //void AddRandomDisplacement();
    //void Integrate();
    //void UpdatePosition();
    Site * GetSite(int i);
    Bond * GetNeighborBond(int i);
    directed_bond GetNeighborDirectedBond(int i);
};

#endif // _SIMCORE_BOND_H_
