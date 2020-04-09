#ifndef _SIMCORE_FILAMENT_GIANT_NUMBER_FLUCTUATION_ANALYSIS_H_
#define _SIMCORE_FILAMENT_GIANT_NUMBER_FLUCTUATION_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"

class GiantNumberFluctuationAnalysis
    : public Analysis<Filament, species_id::filament> {
protected:
  double gnf_box_radius_;
  double gnf_box_increment_;
  int *box_n_;
  RNG *rng_;

  void InitOutput() {
    SetAnalysisName("gnf");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "length diameter bond_length persistence_length umax driving "
               "packing_fraction nsteps nspec delta\n";
    auto it = members_->begin();
    double l = it->GetLength();
    double d = it->GetDiameter();
    double cl = it->GetBondLength();
    double pl = it->GetPersistenceLength();
    double dr = it->GetDriving();
    double nspec = sparams_->n_spec;
    double pf = sparams_->packing_fraction;
    output_ << l << " " << d << " " << cl << " " << pl << " "
            << params_->soft_potential_mag << " " << dr << " " << pf << " "
            << params_->n_steps << " " << nspec << " " << params_->delta
            << "\n";
    int n_boxes = sparams_->number_fluctuation_boxes;
    gnf_box_radius_ =
        0.125 * l / params_->system_radius; // Initial box radius is L/4
    // gnf_box_increment_ = (0.5 - gnf_box_radius_) / (n_boxes - 1);
    for (int j = 0; j < n_boxes; ++j) {
      output_ << "box" << j << " ";
    }
    output_ << "\n";
    double box_rad = gnf_box_radius_;
    for (int j = 0; j < n_boxes; ++j) {
      output_ << box_rad * params_->system_radius << " ";
      box_rad += box_rad;
    }
    output_ << "\n";
    output_ << "time";
    for (int i = 0; i < sparams_->number_fluctuation_centers; ++i) {
      for (int j = 0; j < n_boxes; ++j) {
        output_ << " box" << j;
      }
    }
    output_ << "\n";
    box_n_ = new int[n_boxes];
    rng_ = new RNG(params_->seed);
  }
  void RunAnalysis() {
    int n_boxes = sparams_->number_fluctuation_boxes;
    int n_centers = sparams_->number_fluctuation_centers;
    double center[3] = {0, 0, 0};
    output_ << time_ << " ";
    for (int i_center = 0; i_center < n_centers; ++i_center) {
      for (int i = 0; i < params_->n_dim; ++i) {
        center[i] = rng_->RandomUniform() - 0.5;
      }
      std::fill(box_n_, box_n_ + n_boxes, 0);
      for (auto it = members_->begin(); it != members_->end(); ++it) {
        // Average scaled position of filament
        double asp[3] = {0, 0, 0};
        it->GetAvgScaledPosition(asp);
        bool in_box = false;
        double box_radius = gnf_box_radius_;
        for (int i_box = 0; i_box < n_boxes; ++i_box) {
          if (in_box) {
            box_n_[i_box]++;
            continue;
          } else {
            in_box = true;
          }
          for (int i = 0; i < params_->n_dim; ++i) {
            if (center[i] < 0) {
              if (center[i] - box_radius < -0.5) {
                in_box = in_box && ((asp[i] < center[i] + box_radius) ||
                                    (asp[i] > center[i] - box_radius + 1));
              } else {
                in_box = in_box && ((asp[i] < center[i] + box_radius) &&
                                    (asp[i] > center[i] - box_radius));
              }
            } else {
              if (center[i] + box_radius > 0.5) {
                in_box = in_box && ((asp[i] > center[i] - box_radius) ||
                                    (asp[i] < center[i] + box_radius - 1));
              } else {
                in_box = in_box && ((asp[i] > center[i] - box_radius) &&
                                    (asp[i] < center[i] + box_radius));
              }
            }
          }
          if (in_box) {
            box_n_[i_box]++;
          }
          box_radius += box_radius;
        }
      }
      for (int i_box = 0; i_box < n_boxes; ++i_box) {
        output_ << " " << box_n_[i_box];
      }
    }
    output_ << "\n";
  }

  void EndAnalysis() { 
    delete[] box_n_; 
    delete rng_;
  }
};

#endif // _SIMCORE_FILAMENT_GIANT_NUMBER_FLUCTUATION_ANALYSIS_H_
