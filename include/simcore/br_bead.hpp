#ifndef _SIMCORE_BR_BEAD_H_
#define _SIMCORE_BR_BEAD_H_

#include "object.hpp"

class BrBead : public Object {
protected:
  br_bead_parameters *sparams_;
  bool zero_temperature_ = false;
  double gamma_trans_ = 0;
  double gamma_rot_ = 0;
  double diffusion_ = 0;
  double diffusion_rot_ = 0;
  double driving_factor_ = 0;
  double driving_torque_ = 0;
  int chiral_handedness_ = 0;
  double noise_factor_ = 1;
  bool alignment_interaction_ = false;
  double alignment_torque_ = 0;
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertBrBead();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();

public:
  BrBead(unsigned long seed);
  void Init(br_bead_parameters *sparams);
  void UpdatePosition();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void ZeroForce();
  virtual void ReadSpec(std::fstream &ip);
  virtual void WriteSpec(std::fstream &op);
};
#endif // _SIMCORE_BR_BEAD_H_
