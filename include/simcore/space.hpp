#ifndef _SIMCORE_SPACE_PROPERTIES_H_
#define _SIMCORE_SPACE_PROPERTIES_H_

#include "auxiliary.hpp"
#include "macros.hpp"
#include <math.h>

class Space {
private:
  // geometric data
  int n_dim_ = -1;
  int n_periodic_ = -1;
  boundary_type boundary_;
  std::string boundary_type_string_;
  double radius_ = 0;
  double volume_ = 0;
  double unit_cell_volume_ = 0;
  double a_[9];             // direct lattice vector
  double b_[9];             // reciprocal lattice vector
  double a_perp_[3];        // perp dist between opposite unit cell faces
  double unit_cell_[9];     // unit cell matrix
  double unit_cell_inv_[9]; // inverse unit cell matrix

  // bud data
  double bud_radius_ = 0;  // radius of budding cell
  double bud_height_ = 0;  // center-center dist of mother and daughter cells
  double neck_height_ = 0; // height of mother-daughter cell intersection
  double neck_radius_ = 0; // radius of opening between mother-daughter cells
  double v_ratio_ = 0;     // percentage daughter cell volume to total volume

  // statistical data
  bool constant_pressure_ = false;
  bool constant_volume_ = false;
  bool update_ = false;
  double pressure_ = 0;
  double pressure_tensor_[9];
  double target_pressure_ = 0;
  double target_radius_ = 0;
  double delta_ = 0;
  double compressibility_ = 0; // set to unity (see Berendsen et al. 1984)
  double pressure_time_ = 0;   // time to reach target pressure
  double prev_unit_cell_[9];   // previous unit cell
  double mu_[9];               // scaling matrix for constant pressure

  system_parameters *params_;
  space_struct s_struct;
  void InitUnitCell();
  void InitSpaceStruct();
  void CalculateVolume();
  void CalculateRadius();
  void UpdateVolume();
  void UpdateSpaceStruct();
  void CalculateScalingMatrix();
  void CalculateUnitCellQuantities();
  void UpdateUnitCell();

public:
  Space();
  void Init(system_parameters *params);
  void UpdateSpace();
  void ConstantPressure();
  void ConstantVolume();
  space_struct *GetStruct();
  bool GetUpdate() { return update_; }
};

#endif // _SIMCORE_SPACE_PROPERTIES_H_
