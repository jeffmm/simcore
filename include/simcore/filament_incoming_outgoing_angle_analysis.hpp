#ifndef _SIMCORE_FILAMENT_INCOMING_OUTGOING_ANGLE_ANALYSIS_H_
#define _SIMCORE_FILAMENT_INCOMING_OUTGOING_ANGLE_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class IncomingOutgoingAngleAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  void InitOutput() {
    SetAnalysisName("in_out_angle");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "perlen: " << sparams_->perlen_ratio << "\n";
    output_ << "umax: " << params_->soft_potential_mag << "\n";
    output_ << "angle_in: ";
    if (n_members_ < 2) {
      Logger::Error(
          "Unexpected filament number in IncomingOutgoingAngleAnalysis");
    }
    double u1[3] = {0, 0, 0};
    double u2[3] = {0, 0, 0};
    members_->at(0).GetAvgOrientation(u1);
    members_->at(1).GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim, u1, u2);
    output_ << acos(dp) << "\n";
  }
  void RunAnalysis() {}

  void EndAnalysis() {
    output_ << "angle_out: ";
    double u1[3] = {0, 0, 0};
    double u2[3] = {0, 0, 0};
    members_->at(0).GetAvgOrientation(u1);
    members_->at(1).GetAvgOrientation(u2);
    double dp = dot_product(params_->n_dim, u1, u2);
    output_ << acos(dp) << "\n";
  }
};

#endif // _SIMCORE_FILAMENT_INCOMING_OUTGOING_ANGLE_ANALYSIS_H_
