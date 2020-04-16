#ifndef _SIMCORE_FILAMENT_MSD_ANALYSIS_H_
#define _SIMCORE_FILAMENT_MSD_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class MSDAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  bool use_center_of_curvature_ = false;
  void InitOutput() {
    SetAnalysisName("msd");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "nsteps nspec delta n_filaments length diameter bond_length "
               "persistence_length driving curvature\n";
    auto it = members_->begin();
    double nspec = sparams_->n_spec;
    double l = sparams_->length;
    double lp = sparams_->persistence_length;
    if (sparams_->perlen_ratio > 0) {
      lp = sparams_->perlen_ratio * l;
    }
    double dr = sparams_->driving_factor;
    if (sparams_->peclet_number > 0) {
      dr = sparams_->peclet_number / SQR(l);
    }
    if (sparams_->polydispersity_flag || sparams_->dynamic_instability_flag) {
      Logger::Warning("Filament MSD analysis assumes monodisperse filaments and"
                      " no dynamic instability. Results may be nonsense.");
    }
    double ic = sparams_->intrinsic_curvature;
    if (sparams_->radius_of_curvature > 0) {
      ic = 1.0 / sparams_->radius_of_curvature;
    }
    if (ic > 0) {
      use_center_of_curvature_ = true;
    }
    output_ << params_->n_steps << " " << sparams_->n_spec << " "
            << params_->delta << " " << n_members_ << " " << l << " "
            << sparams_->diameter << " " << members_->begin()->GetBondLength()
            << " " << lp << " " << dr << " " << ic << "\n";
    output_ << "time";
    for (int i = 0; i < n_members_; ++i) {
      output_ << " pos_x pos_y pos_z u_x u_y u_z";
    }
    output_ << "\n";
  }
  void RunAnalysis() {
    double avg_pos[3] = {0};
    double avg_u[3] = {0};
    double avg_spos[3] = {0};
    double cc[3] = {0};
    output_ << time_;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      it->GetAvgPosition(avg_pos);
      if (use_center_of_curvature_) {
        // Use vector from avg pos to center of curvature for orientation vector
        it->GetAvgScaledPosition(avg_spos);
        it->GetCenterOfCurvature(cc);
        for (int i = 0; i < params_->n_periodic; ++i) {
          avg_u[i] = cc[i] - avg_spos[i];
          avg_u[i] -= NINT(avg_u[i]);
        }
        normalize_vector(avg_u, params_->n_dim);
      } else {
        it->GetAvgOrientation(avg_u);
      }
      for (int i = 0; i < 3; ++i) {
        output_ << " " << avg_pos[i];
      }
      for (int i = 0; i < 3; ++i) {
        output_ << " " << avg_u[i];
      }
    }
    output_ << "\n";
    output_ << std::flush;
  }

  void EndAnalysis() {
    /* output a handedness file if we have curved filaments */
    if (use_center_of_curvature_) {
      std::string file_name = params_->run_name + "_" + sid_._to_string() +
                              "_" + sparams_->name + ".handedness.analysis";

      Logger::Debug("Writing analysis output file %s", file_name.c_str());
      std::ofstream hand_out(file_name);
      hand_out << "n_filaments intrinsic_curvature\n";
      hand_out << n_members_ << " " << sparams_->intrinsic_curvature << "\n";
      hand_out << "handedness\n";
      for (auto it=members_->begin(); it!=members_->end(); ++it) {
        hand_out << it->GetHandedness() << " ";
      }
      hand_out << "\n";
      hand_out.close();
    }
  }
};

#endif // _SIMCORE_FILAMENT_MSD_ANALYSIS_H_
