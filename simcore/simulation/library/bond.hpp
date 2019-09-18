#ifndef _SIMCORE_BOND_H_
#define _SIMCORE_BOND_H_

#include "site.hpp"

// Bonds, ie graph edges
class Bond : public Object {
protected:
  Site *sites_[2];
  int bond_number_ = 0;
  double dr_[3], equil_length_, k_spring_, body_frame_[6], orientation_0_[3];
  Object *mesh_ptr_;

public:
  Bond();
  void Init(Site *s1, Site *s2);
  void ReInit();
  void Report();
  void ReportSites();
  void Draw(std::vector<graph_struct *> *graph_array);
  int const GetBondNumber();
  void SetBondNumber(int bnum);
  void ZeroOrientationCorrelation();
  double const GetOrientationCorrelation();
  // void GetBodyFrame();
  // void UpdateOrientation();
  // void AddRandomDisplacement();
  // void Integrate();
  // void UpdatePosition();
  Site *GetSite(int i);
  Bond *GetNeighborBond(int i);
  directed_bond GetNeighborDirectedBond(int i);
  virtual bool HasNeighbor(int other_oid);
  void SetMeshPtr(Object * obj_ptr);
  Object * GetMeshPtr() { return mesh_ptr_; }
};

#endif // _SIMCORE_BOND_H_
