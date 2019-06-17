#ifndef _SIMCORE_ANCHOR_H_
#define _SIMCORE_ANCHOR_H_

#include "bond.h"

class Anchor : public Object {
  private:
    bool bound_,
         walker_,
         diffuse_,
         active_;

    int step_direction_;

    double bond_length_,
           bond_lambda_,
           mesh_lambda_,
           velocity_,
           max_velocity_,
           diffusion_,
           f_spring_max_;

    Bond * bond_;

  public:
    Anchor();
    void Init();
    void UpdatePosition();
    void Diffuse();
    //void DiffuseBound();
    void Activate();
    void Deactivate();
    void CheckNearBoundary();
    void CheckNearBuddingBoundary();
    void AnchorBoundary(double * anchor_point);
    void DetachBoundary();
    void ApplyAnchorForces();
    void AttachToBond(directed_bond, double lambda, double mesh_lambda);
    bool SwitchBonds(bool next_bond, double lambda);
    void UpdateAnchorPosition() {}
    void SetDiffusion();
    void SetWalker(int dir,double walk_v);
    void Walk();
    void AttachBondRandom(Bond * b, double mesh_lambda);
    double const GetMeshLambda();
    void Clear();
};

#endif
