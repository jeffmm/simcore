#ifndef _SIMCORE_SPHEROCYLINDER_H_
#define _SIMCORE_SPHEROCYLINDER_H_

#include "br_rod.hpp"
class Spherocylinder : public BrRod {
 protected:
  spherocylinder_parameters *sparams_;
  void ApplyForcesTorques();
  void InsertSpherocylinder();
  void SetParameters();

 public:
  Spherocylinder(unsigned long seed);
  void Init(spherocylinder_parameters *sparams);
  void UpdatePosition();
};

#endif  // _SIMCORE_SPHEROCYLINDER_H_
