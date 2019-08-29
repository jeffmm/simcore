#ifndef _SIMCORE_SPACE_PROPERTIES_H_
#define _SIMCORE_SPACE_PROPERTIES_H_

#include <math.h>
#include "auxiliary.hpp"
#include "macros.hpp"

class Space {
 private:
  // geometric data
  int n_dim_, n_periodic_;
  boundary_type boundary_;
  std::string boundary_type_string_;
  double radius_, volume_, unit_cell_volume_,
      a_[9],       // direct lattice vector
      b_[9],       // reciprocal lattice vector
      a_perp_[3],  // perp dist between opposite unit cell faces
      unit_cell_[9],
      unit_cell_inv_[9];  // inverse unit cell matrix

  // bud data
  double bud_radius_,  // radius of budding cell
      bud_height_,     // center-center separation of parent and daughter cells
      neck_height_,    // height of plane at parent-daughter cell intersection
      neck_radius_,    // radius of circle defined at parent-daughter cell
                       // intersection
      v_ratio_;        // volume ratio of daughter cell to total system volume

  // statistical data
  bool constant_pressure_, constant_volume_, update_;
  double pressure_, pressure_tensor_[9], target_pressure_, target_radius_,
      delta_,
      compressibility_,    // set to unity (see Berendsen et al. 1984)
      pressure_time_,      // time to reach target pressure
      prev_unit_cell_[9],  // previous unit cell
      mu_[9];              // scaling matrix for constant pressure

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

#endif  // _SIMCORE_SPACE_PROPERTIES_H_
