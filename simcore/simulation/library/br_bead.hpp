#ifndef _SIMCORE_BR_BEAD_H_
#define _SIMCORE_BR_BEAD_H_

#include "object.hpp"

class BrBead : public Object {
 protected:
  bool stoch_flag_;
  double gamma_trans_, gamma_rot_, diffusion_, driving_factor_;
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertBrBead();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();

 public:
  BrBead();
  void Init();
  void UpdatePosition();
  virtual void GetInteractors(std::vector<Object *> *ix);
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> *graph_array);
  virtual void ZeroForce();
};
#endif  // _SIMCORE_BR_BEAD_H_
