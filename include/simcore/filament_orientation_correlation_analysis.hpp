#ifndef _SIMCORE_FILAMENT_ORIENTATION_CORRELATION_ANALYSIS_H_
#define _SIMCORE_FILAMENT_ORIENTATION_CORRELATION_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class OrientationCorrelationAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  int n_iterations_;
  void InitOutput() {
    SetAnalysisName("orientation_correlation");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    n_iterations_ = sparams_->orientation_corr_n_steps;
    output_ << "n_filaments n_bonds n_avg_steps\n";
    output_ << n_members_ << members_->begin()->GetNBonds() << n_iterations_
            << "\n";
    output_ << "time orientation_corr_avg orientation_corr_sem\n";
  }
  void RunAnalysis() {
    if (iteration_ % n_iterations_ == 0) {
      for (auto it = members_->begin(); it != members_->end(); ++it) {
        it->ZeroOrientationCorrelations();
      }
      output_ << "0 1 0\n";
      return;
    }
    double avg = 0;
    double sem = 0;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      std::pair<double, double> avg_var = it->GetAvgOrientationCorrelation();
      avg += avg_var.first;
      sem += avg_var.second;
    }
    avg /= n_members_;
    sem = sem / n_members_; // This calculates the pooled variance
    sem = sqrt(sem / (n_members_ *
                      members_->begin()->GetNBonds())); // """ standard error
    output_ << 0.5 * (iteration_ % n_iterations_) * params_->delta << " " << avg
            << " " << sem << "\n";
  }

  void EndAnalysis() {}
};

#endif // _SIMCORE_FILAMENT_ORIENTATION_CORRELATION_ANALYSIS_H_
