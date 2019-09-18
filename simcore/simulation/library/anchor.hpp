#ifndef _SIMCORE_ANCHOR_H_
#define _SIMCORE_ANCHOR_H_

#include "mesh.hpp"

class Anchor : public Object {
private:
  bool bound_;
  bool walker_;
  bool diffuse_;
  bool active_;
  bool end_pausing_;

  int step_direction_;

  double bond_length_;
  double bond_lambda_;
  double mesh_lambda_;
  double velocity_;
  double max_velocity_, diffusion_;
  double f_spring_max_;
  double k_off_;
  double f_stall_;
  double force_dep_vel_flag_;

  Bond *bond_;
  Mesh *mesh_;

  int bond_oid_;
  int mesh_n_bonds_;

public:
  Anchor();
  void Init();
  bool IsBound();
  void UpdatePosition();
  void Diffuse();
  // void DiffuseBound();
  void Activate();
  void Deactivate();
  // void CheckNearBoundary();
  // void CheckNearBuddingBoundary();
  // void AnchorBoundary(double * anchor_point);
  // void DetachBoundary();
  void ApplyAnchorForces();
  void AttachToBond(directed_bond, double lambda, double mesh_lambda);
  bool SwitchBonds(bool next_bond, double lambda);
  void UpdateAnchorPosition() {}
  void SetDiffusion();
  void SetWalker(int dir, double walk_v);
  void Walk();
  void AttachObjRandom(Object *o);
  void AttachObjLambda(Object *o, double lambda);
  double const GetMeshLambda();
  double const GetBondLambda();
  void SetBondLambda(double l);
  void SetBound();
  void Clear();
  int const GetBoundOID();
  void Draw(std::vector<graph_struct *> *graph_array);
};

#endif
