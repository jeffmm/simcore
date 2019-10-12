#ifndef _SIMCORE_ANCHOR_H_
#define _SIMCORE_ANCHOR_H_

#include "mesh.hpp"
#include "neighbor_list.hpp"

/* Class for bound crosslink heads (called anchors). Tracks and updates its
   absolute position in space and relative position to its bound object. */
class Anchor : public Object {
private:
  bool bound_;
  bool walker_;
  bool diffuse_;
  bool static_flag_;
  bool active_;
  bool end_pausing_;
  crosslink_parameters *sparams_;
  int step_direction_;

  double bond_length_;
  double bond_lambda_;
  double mesh_length_;
  double mesh_lambda_;
  double velocity_;
  double max_velocity_, diffusion_;
  double k_off_;
  double f_stall_;
  double force_dep_vel_flag_;

  NeighborList neighbors_;

  Bond *bond_;
  Mesh *mesh_;

  int mesh_n_bonds_;

  void UpdateAnchorPositionToBond();
  void Diffuse();
  void Walk();
  bool CheckMesh();
  bool CalcBondLambda();

public:
  Anchor();
  void Init(crosslink_parameters *sparams);
  bool IsBound();
  void UpdatePosition();
  void Activate();
  void Deactivate();
  void ApplyAnchorForces();
  void UpdateAnchorPositionToMesh();
  void SetDiffusion();
  void SetWalker(int dir, double walk_v);
  void AttachObjRandom(Object *o);
  void AttachObjLambda(Object *o, double lambda);
  void AttachObjMeshLambda(Object *o, double mesh_lambda);
  double const GetMeshLambda();
  double const GetBondLambda();
  void SetBondLambda(double l);
  void SetMeshLambda(double ml);
  void SetBound();
  void Unbind();
  int const GetBoundOID();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighbor(Object *neighbor);
  void ClearNeighbors();
  const Object *const *GetNeighborListMem();
  Object *GetNeighbor(int i_neighbor);
  const int GetNNeighbors() const;
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void BindToPosition(double *pos);
  void SetStatic(bool static_flag);
};

#endif
