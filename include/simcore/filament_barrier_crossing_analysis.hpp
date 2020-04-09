#ifndef _SIMCORE_FILAMENT_CROSSING_ANALYSIS_H_
#define _SIMCORE_FILAMENT_CROSSING_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class BarrierCrossingAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  void InitOutput() {
    SetAnalysisName("barrier_crossing");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_
        << "barrier_resistance persistence_length_ratio driving n_filaments\n";
    double lp_ratio = sparams_->perlen_ratio;
    if (lp_ratio < 0) {
      lp_ratio =
          members_->at(0).GetPersistenceLength() / members_->at(0).GetLength();
    }
    if (params_->potential.compare("wca") == 0) {
      output_ << "INF ";
    } else {
      output_ << params_->soft_potential_mag << " ";
    }
    output_ << lp_ratio << " " << sparams_->driving_factor << " " << n_members_
            << "\n";

    double avg_u[3] = {0, 0, 0};
    double avg_pos[3] = {0, 0, 0};
    output_ << "initial_angles_rel_to_barrier\n";
    bool do_warning = false;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      it->GetAvgOrientation(avg_u);
      it->GetAvgPosition(avg_pos);
      output_ << acos(avg_u[1]) << " ";
      if (avg_pos[0] > 0) {
        do_warning = true;
      }
    }
    output_ << "\n";
    output_ << "warning: filaments initialized in an unexpected configuration\n";

    if (do_warning) {
      Logger::Warning("CrossingAnalysis assumes filaments are initialized on"
                      " the left-hand side of the simulation box. The results"
                      " may be nonsense.");
    }
  }
  void RunAnalysis() {}

  void EndAnalysis() {
    double avg_pos[3] = {0, 0, 0};
    output_ << "crossed_barrier\n";
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      it->GetAvgPosition(avg_pos);
      output_ << (avg_pos[0] > 2 ? 1 : 0) << " ";
    }
    output_ << "\n";
  }
};

#endif // _SIMCORE_FILAMENT_CROSSING_ANALYSIS_H_
