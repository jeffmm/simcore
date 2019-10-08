#ifndef _SIMCORE_SPINDLE_H_
#define _SIMCORE_SPINDLE_H_

#include "anchor.hpp"
#include "filament.hpp"
#include "spherocylinder.hpp"

class Spindle : public Spherocylinder {
protected:
  spindle_parameters *sparams_;
  bool alignment_potential_, fixed_spacing_;
  int n_filaments_bud_, n_filaments_mother_;
  double k_spring_, k_align_, spring_length_, anchor_distance_, gamma_trans_,
      gamma_rot_, diffusion_, spb_diameter_;
  std::vector<Filament> filaments_;
  std::vector<Anchor> anchors_;
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertSpindle();
  void GenerateAnchorSites();
  void Integrate();
  void ResetAnchorPositions();
  bool InsertFilament(int i);
  void SetParameters();

public:
  Spindle();
  void Init(spindle_parameters *sparams);
  void UpdatePosition() {}
  void UpdatePosition(bool midstep);
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void ZeroForce();
};

#endif // _SIMCORE_SPINDLE_H_
