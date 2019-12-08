#ifndef _SIMCORE_BR_ROD_H_
#define _SIMCORE_BR_ROD_H_

#include "object.hpp"

class BrRod : public Object {
 protected:
  double gamma_par_ = 0;
  double gamma_perp_ = 0;
  double gamma_rot_ = 0;
  double diffusion_par_ = 0;
  double diffusion_perp_ = 0;
  double diffusion_rot_ = 0;
  double body_frame_[6];
  void InsertRod(std::string insertion_type, double buffer = -1);
  void SetDiffusion();
  void GetBodyFrame();
  void AddRandomDisplacement();
  void AddRandomReorientation();
  void Integrate();

 public:
  BrRod(unsigned long seed);
};

#endif  // _SIMCORE_BR_ROD_H_
