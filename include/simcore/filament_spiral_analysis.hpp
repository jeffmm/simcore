#ifndef _SIMCORE_FILAMENT_SPIRAL_ANALYSIS_H_
#define _SIMCORE_FILAMENT_SPIRAL_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class SpiralAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  void InitOutput() {
    SetAnalysisName("spiral");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "spiral_analysis_file\n";
    output_ << "n_filaments nsteps nspec delta";
    for (int i = 0; i < n_members_; ++i) {
      output_ << " length diameter bond_length persistence_length driving";
    }
    output_ << "\n";
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      double l = it->GetLength();
      double d = it->GetDiameter();
      double cl = it->GetBondLength();
      double pl = it->GetPersistenceLength();
      double dr = it->GetDriving();
      double nspec = sparams_->n_spec;
      output_ << l << " " << d << " " << cl << " " << pl << " " << dr << " "
              << params_->n_steps << " " << nspec << " " << params_->delta
              << "\n";
    }
    output_ << "time";
    for (int i = 0; i < n_members_; ++i) {
      output_ << " angle_sum E_bend tip_z_proj spiral_number head_pos_x "
                 "head_pos_y tail_pos_x tail_pos_y";
    }
    output_ << "\n";
  }
  void RunAnalysis() {
    // Analysis here
    output_ << time_;
    for (auto it=members_->begin(); it!=members_->end(); ++it) {
      double e_bend = 0;
      double tot_angle = 0;
      double length = it->GetLength();
      double plength = it->GetPersistenceLength();
      double clength = it->GetBondLength();
      double e_zero = length * plength / (clength * clength);
      it->CalculateSpiralNumber();
      double spiral_number = it->GetSpiralNumber();
      std::vector<double> const *const thetas = it->GetThetas();
      for (int i = 0; i < thetas->size(); ++i) {
        tot_angle += acos((*thetas)[i]);
        e_bend += (*thetas)[i];
      }
      // record energy relative to the bending "zero energy" (straight rod)
      e_bend = e_zero - e_bend * plength / clength;
      double tip_z = it->GetTipZ();
      double const *const head_pos = it->GetHeadPosition();
      double const *const tail_pos = it->GetTailPosition();
      output_ << " " << tot_angle << " " << e_bend << " " << tip_z << " "
              << spiral_number << " " << head_pos[0] << " " << head_pos[1] << " "
              << tail_pos[0] << " " << tail_pos[1];
    }
    output_ << "\n";
  }

  void EndAnalysis() {}
};

#endif // _SIMCORE_FILAMENT_SPIRAL_ANALYSIS_H_
