#ifndef _SIMCORE_FILAMENT_SPECIES_H_
#define _SIMCORE_FILAMENT_SPECIES_H_

#include "filament.hpp"
#include "species.hpp"
#include "filament_curvature_cluster_analysis.hpp"
#include "filament_spiral_analysis.hpp"
#include "filament_angle_distribution_analysis.hpp"
#include "filament_flocking_analysis.hpp"
#include "filament_orientation_correlation_analysis.hpp"
#include "filament_polar_order_analysis.hpp"
#include "filament_global_order_analysis.hpp"
#include "filament_end_to_end_fluctuation_analysis.hpp"
#include "filament_giant_number_fluctuation_analysis.hpp"
#include "filament_barrier_crossing_analysis.hpp"
#include "filament_incoming_outgoing_angle_analysis.hpp"
#include "filament_mean_squared_displacement_analysis.hpp"

typedef Analysis<Filament, species_id::filament> FilamentAnalysis;

class FilamentSpecies : public Species<Filament, species_id::filament> {
 protected:
  bool midstep_ = true;
  // Analysis structures
  double fill_volume_;
  double packing_fraction_;
  std::fstream error_file_;
  void LoadAnalysis();
  void InitErrorAnalysis();
  void RunErrorAnalysis();

 public:
  FilamentSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void ResetPreviousPositions();
  void AddMember();

  void Reserve();
  void UpdatePositions();
  void CleanUp();
  virtual const double GetSpecLength() const;
  // Redundant for filaments.
  virtual void CenteredOrientedArrangement() {}

};

#endif
