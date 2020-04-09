#ifndef _SIMCORE_FILAMENT_GLOBAL_ORDER_ANALYSIS_H_
#define _SIMCORE_FILAMENT_GLOBAL_ORDER_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class GlobalOrderAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  double *nematic_order_tensor_;
  double *polar_order_vector_;

  void InitOutput() {
    SetAnalysisName("global_order");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "global_order_analysis_file\n";
    output_
        << "time polar_order_x polar_order_y polar_order_z nematic_order_xx "
           "nematic_order_xy nematic_order_xz nematic_order_yx nematic_order_yy "
           "nematic_order_yz nematic_order_zx nematic_order_zy nematic_order_zz "
           "spiral_order signed_spiral_order n_spooling "
           "avg_spool_spiral_number\n";
    nematic_order_tensor_ = new double[9];
    polar_order_vector_ = new double[3];
    std::fill(nematic_order_tensor_, nematic_order_tensor_ + 9, 0.0);
    std::fill(polar_order_vector_, polar_order_vector_ + 3, 0.0);
  }
  void RunAnalysis() {
    double sn_tot = 0.0;
    double sn_mag = 0.0;
    double sn_spools = 0;
    int n_spooling = 0;
    double sn;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      it->CalculateSpiralNumber();
      sn = it->GetSpiralNumber();
      it->GetPolarOrder(polar_order_vector_);
      it->GetNematicOrder(nematic_order_tensor_);
      sn_mag += ABS(sn);
      sn_tot += sn;
      if (sn > 0.7) {
        sn_spools += sn;
        n_spooling++;
      }
    }
    sn_mag /= n_members_;
    sn_tot /= n_members_;
    if (n_spooling > 0) {
      sn_spools /= n_spooling;
    } else {
      sn_spools = 0;
    }
    for (int i = 0; i < 3; ++i) {
      polar_order_vector_[i] /= n_members_;
    }
    for (int i = 0; i < 9; ++i) {
      nematic_order_tensor_[i] /= n_members_;
    }
    if (output_.is_open()) {
      output_ << time_ << " " << polar_order_vector_[0] << " "
                         << polar_order_vector_[1] << " "
                         << polar_order_vector_[2] << " "
                         << nematic_order_tensor_[0] << " "
                         << nematic_order_tensor_[1] << " "
                         << nematic_order_tensor_[2] << " "
                         << nematic_order_tensor_[3] << " "
                         << nematic_order_tensor_[4] << " "
                         << nematic_order_tensor_[5] << " "
                         << nematic_order_tensor_[6] << " "
                         << nematic_order_tensor_[7] << " "
                         << nematic_order_tensor_[8] << " " << sn_mag << " "
                         << sn_tot << " " << n_spooling << " " << sn_spools
                         << "\n";
    } else {
      Logger::Error("Problem opening output file in GlobalOrderAnalysis!");
    }
  }

  void EndAnalysis() {
    delete[] nematic_order_tensor_;
    delete[] polar_order_vector_;
  }
};

#endif // _SIMCORE_FILAMENT_GLOBAL_ORDER_ANALYSIS_H_
