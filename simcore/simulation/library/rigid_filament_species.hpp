#ifndef _SIMCORE_RIGID_FILAMENT_SPECIES_H_
#define _SIMCORE_RIGID_FILAMENT_SPECIES_H_

#include "rigid_filament.hpp"
#include "species.hpp"

class RigidFilamentSpecies
    : public Species<RigidFilament, species_id::rigid_filament> {
 protected:
  bool midstep_ = true;
  // Analysis structures
  double fill_volume_;
  double packing_fraction_;
  double e_bend_;
  double tot_angle_;
  double mse2e_;
  double mse2e2_;
  double polar_bin_width_;
  double contact_bin_width_;
  double contact_cut_;
  double po_avg_;
  double cn_avg_;
  double *nematic_order_tensor_;
  double *polar_order_vector_;
  int **theta_histogram_;
  int *polar_order_histogram_;
  int time_;
  int n_bins_;
  int n_bins_1d_;
  int n_samples_;
  int orientation_corr_n_steps_;
  std::fstream spiral_file_;
  std::fstream flock_file_;
  std::fstream theta_file_;
  std::fstream mse2e_file_;
  std::fstream global_order_file_;
  std::fstream local_order_file_;
  std::fstream polar_order_file_;
  std::fstream orientation_corr_file_;
  std::fstream crossing_file_;
  std::fstream polar_order_avg_file_;
  std::fstream in_out_file_;

 public:
  RigidFilamentSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();

  void AddMember();

  void Reserve();
  void UpdatePositions();
  // Redundant for filaments.
  virtual void CenteredOrientedArrangement() {}

  virtual const double GetSpecLength() const;

  void InitAnalysis();
  void RunAnalysis();
  void FinalizeAnalysis();

  void InitSpiralAnalysis();
  void RunSpiralAnalysis();
  // void FinalizeSpiralAnalysis();

  void InitCrossingAnalysis();
  void RunCrossingAnalysis();
  void FinalizeCrossingAnalysis();

  void InitThetaAnalysis();
  void RunThetaAnalysis();
  void FinalizeThetaAnalysis();

  void InitMse2eAnalysis();
  void RunMse2eAnalysis();
  void FinalizeMse2eAnalysis();

  void InitGlobalOrderAnalysis();
  void RunGlobalOrderAnalysis();
  void FinalizeGlobalOrderAnalysis();

  void InitPolarOrderAnalysis();
  void RunPolarOrderAnalysis();
  void FinalizePolarOrderAnalysis();

  void InitLocalOrderAnalysis();
  void RunLocalOrderAnalysis();
  void WriteLocalOrderData();
  void FinalizeLocalOrderAnalysis();

  void InitOrientationCorrelationAnalysis();
  void RunOrientationCorrelationAnalysis();
  void FinalizeOrientationCorrelationAnalysis();

  void InitFlockingAnalysis();
  void RunFlockingAnalysis();
  void FinalizeFlockingAnalysis();
};

#endif
