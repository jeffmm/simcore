#ifndef _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_
#define _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class CurvatureClusterAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  double **curvature_centers_;
  double *curvature_radii_;

  void InitOutput() {
    SetAnalysisName("curve_cluster");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    curvature_centers_ = new double *[n_members_];
    curvature_radii_ = new double[n_members_];
    std::fill(curvature_radii_, curvature_radii_ + 3, 0.0);
    for (int i = 0; i < n_members_; ++i) {
      curvature_centers_[i] = new double[3];
      std::fill(curvature_centers_[i], curvature_centers_[i] + 3, 0.0);
    }
    output_ << "Filler text for now\n";
  }
  void RunAnalysis() {
    // Analysis here
  }

  void EndAnalysis() {
    for (int i = 0; i < n_members_; ++i) {
      delete [] curvature_centers_[i];
    }
    delete [] curvature_centers_;
    delete [] curvature_radii_;
  }
};

#endif // _SIMCORE_FILAMENT_CURVATURE_CLUSTER_ANALYSIS_H_
