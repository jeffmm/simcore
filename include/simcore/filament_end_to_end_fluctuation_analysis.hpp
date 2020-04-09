#ifndef _SIMCORE_FILAMENT_END_TO_END_FLUCTUATION_ANALYSIS_H_
#define _SIMCORE_FILAMENT_END_TO_END_FLUCTUATION_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class EndToEndFluctuationAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  double mse2e_;
  double mse2e2_;
  int n_samples_;

  void InitOutput() {
    SetAnalysisName("mse2e");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "mean_square_end_to_end_analysis\n";
    output_ << "length diameter bond_length persistence_length driving ndim "
                   "nsteps nspec delta theory\n";
    auto it = members_->begin();
    double l = it->GetLength();
    double d = it->GetDiameter();
    double cl = it->GetBondLength();
    double pl = it->GetPersistenceLength();
    double dr = it->GetDriving();
    double nspec = sparams_->n_spec;
    double theory;
    if (params_->n_dim == 2) {
      theory = l * pl * 4.0 - 8.0 * pl * pl * (1 - exp(-0.5 * l / pl));
    } else {
      theory = l * pl * 2.0 - 2.0 * pl * pl * (1 - exp(-l / pl));
    }
    output_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
                << params_->n_dim << " " << params_->n_steps << " " << nspec
                << " " << params_->delta << " " << theory << "\n";
    output_ << "num_filaments_averaged mse2e_mean mse2e_std_err\n";
    mse2e_ = 0.0;
    mse2e2_ = 0.0;
    n_samples_ = 0;
  }
  void RunAnalysis() {
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      double const *const head_pos = it->GetHeadPosition();
      double const *const tail_pos = it->GetTailPosition();
      double mse2e_temp = 0.0;
      for (int i = 0; i < params_->n_dim; ++i) {
        double temp = (head_pos[i] - tail_pos[i]);
        mse2e_temp += temp * temp;
      }
      mse2e_ += mse2e_temp;
      mse2e2_ += mse2e_temp * mse2e_temp;
    }
    n_samples_++;
  }

  void EndAnalysis() {
    output_ << n_members_ << " ";
    mse2e_ /= n_samples_ * n_members_;
    mse2e2_ /= n_samples_ * n_members_;
    output_ << mse2e_ << " ";
    output_ << sqrt((mse2e2_ - mse2e_ * mse2e_) / (n_members_ * n_samples_)) << "\n";
  }
};

#endif // _SIMCORE_FILAMENT_END_TO_END_FLUCTUATION_ANALYSIS_H_
