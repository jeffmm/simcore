#ifndef _SIMCORE_SPHEROCYLINDER_H_
#define _SIMCORE_SPHEROCYLINDER_H_

#include "object.hpp"
class Spherocylinder : public Object {
protected:
  spherocylinder_parameters *sparams_;
  double gamma_par_, gamma_perp_, gamma_rot_, diffusion_par_, diffusion_perp_,
      diffusion_rot_, body_frame_[6];
  bool is_midstep_;
  void ApplyForcesTorques();
  void InsertSpherocylinder();
  void SetDiffusion();
  void GetBodyFrame();
  void AddRandomDisplacement();
  void AddRandomReorientation();
  void Integrate();
  void SetParameters();

public:
  Spherocylinder();
  void Init(spherocylinder_parameters *sparams);
  void UpdatePosition();
};

#endif // _SIMCORE_SPHEROCYLINDER_H_
