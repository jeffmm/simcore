#ifndef _SIMCORE_FILAMENT_ANGLE_DISTRIBUTION_ANALYSIS_H_
#define _SIMCORE_FILAMENT_ANGLE_DISTRIBUTION_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class AngleDistributionAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  int n_bins_;
  int **theta_histogram_;
  void InitOutput() {
    SetAnalysisName("angle_distribution");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    if (params_->interaction_flag) {
      Logger::Warning("Conducting angle distribution analysis on interacting"
                      " filaments");
    }
    output_
        << "length diameter bond_length persistence_length driving n_filaments "
           "n_bonds n_steps n_spec delta n_dim \n";
    double l, cl, pl, dr, d;
    int nbonds;
    int nmembers = members_->size();
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      l = it->GetLength();
      d = it->GetDiameter();
      cl = it->GetBondLength();
      pl = it->GetPersistenceLength();
      dr = it->GetDriving();
      nbonds = it->GetNBonds();
    }
    int nspec = sparams_->n_spec;
    output_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
            << nmembers << " " << nbonds << " " << params_->n_steps << " "
            << nspec << " " << params_->delta << " " << params_->n_dim << "\n";
    output_ << "cos_theta";
    for (int i = 0; i < nbonds - 1; ++i) {
      output_ << " theta_" << i + 1 << i + 2;
    }
    output_ << "\n";
    n_bins_ = 10000;
    theta_histogram_ = new int *[nbonds - 1];
    for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
      theta_histogram_[ibond] = new int[n_bins_];
      for (int ibin = 0; ibin < n_bins_; ++ibin) {
        theta_histogram_[ibond][ibin] = 0;
      }
    }
  }

  void RunAnalysis() {
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      std::vector<double> const *const thetas = it->GetThetas();
      for (int i = 0; i < (it->GetNBonds() - 1); ++i) {
        int bin_number = (int)floor((1 + (*thetas)[i]) * (n_bins_ / 2));
        if (bin_number == n_bins_) {
          bin_number = n_bins_ - 1;
        } else if (bin_number == -1) {
          bin_number = 0;
        } else if (bin_number > n_bins_ && bin_number < 0) {
          Logger::Error("Something went wrong in AngleDistributionAnalysis!");
        }
        theta_histogram_[i][bin_number]++;
      }
    }
  }

  void EndAnalysis() {
    int nbonds = members_->back().GetNBonds();
    for (int i = 0; i < n_bins_; ++i) {
      double axis = (2.0 / n_bins_) * i - 1;
      output_ << " " << axis;
      for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
        output_ << " " << theta_histogram_[ibond][i];
      }
      output_ << "\n";
    }

    for (int ibond = 0; ibond < nbonds - 1; ++ibond) {
      delete[] theta_histogram_[ibond];
    }
    delete[] theta_histogram_;
  }
};

#endif // _SIMCORE_FILAMENT_ANGLE_DISTRIBUTION_ANALYSIS_H_
