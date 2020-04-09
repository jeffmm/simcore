#ifndef _SIMCORE_FILAMENT_FLOCKING_ANALYSIS_H_
#define _SIMCORE_FILAMENT_FLOCKING_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class FlockingAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  int *flock_states_;
  void InitOutput() {
    SetAnalysisName("flock");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "length diameter bond_length persistence_length umax driving "
               "nsteps nspec delta\n";
    auto it = members_->begin();
    double l = it->GetLength();
    double d = it->GetDiameter();
    double cl = it->GetBondLength();
    double pl = it->GetPersistenceLength();
    double dr = it->GetDriving();
    double nspec = sparams_->n_spec;
    output_ << l << " " << d << " " << cl << " " << pl << " "
            << params_->soft_potential_mag << " " << dr << " "
            << params_->n_steps << " " << nspec << " " << params_->delta
            << "\n";
    output_ << "time n_flocking n_exterior n_interior n_joined n_left";
    for (int i = 0; i < n_members_; ++i) {
      output_ << " fil" << i;
    }
    output_ << "\n";
    flock_states_ = new int[n_members_];
    RequireInteractionAnalysis();
  }
  void RunAnalysis() {
    int n_flocking = 0;
    int n_interior = 0;
    int n_exterior = 0;
    int n_joined = 0;
    int n_left = 0;
    int i = 0;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      it->CheckFlocking();
      int flock_type = it->GetFlockType();
      int flock_change_state = it->GetFlockChangeState();
      flock_states_[i++] = flock_type;
      if (flock_type) {
        n_flocking++;
        if (flock_type == 1) {
          n_exterior++;
        } else if (flock_type == 2) {
          n_interior++;
        } else {
          Logger::Warning("Flock type not recognized in "
                          "FlockingAnalysis");
        }
      }
      if (flock_change_state) {
        if (flock_change_state == 1) {
          n_joined++;
        } else if (flock_change_state == 2) {
          n_left++;
        } else {
          Logger::Warning("Flock change state not recognized in "
                          "FlockingAnalysis");
        }
      }
    }
    if (output_.is_open()) {
      output_ << time_ << " " << n_flocking << " " << n_exterior << " "
              << n_interior << " " << n_joined << " " << n_left;
      for (int i = 0; i < n_members_; ++i) {
        output_ << " " << flock_states_[i];
      }
      output_ << "\n";
    } else {
      Logger::Error("Problem opening file in FlockingAnalysis!");
    }
  }

  void EndAnalysis() { delete[] flock_states_; }
};

#endif // _SIMCORE_FILAMENT_FLOCKING_ANALYSIS_H_
