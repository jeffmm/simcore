#ifndef _SIMCORE_FILAMENT_POLAR_ORDER_ANALYSIS_H_
#define _SIMCORE_FILAMENT_POLAR_ORDER_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"
#include "cell_list.hpp"

class PolarOrderAnalysis : public Analysis<Filament, species_id::filament> {
protected:
  int n_bins_;
  int n_bins_1d_;
  int *histogram_;
  double polar_bin_width_;
  double contact_bin_width_;
  double contact_cut_;
  std::fstream output_avg_;

  void InitOutput() {
    SetAnalysisName("polar_order");
    Analysis::InitOutput();
    SetAnalysisName("polar_order_avg");
    std::string file_name = params_->run_name + "_" + sid_._to_string() + "_" +
                            sparams_->name + "." + GetAnalysisName() +
                            ".analysis";
    Logger::Debug("Initializing analysis output file %s", file_name.c_str());
    output_avg_.open(file_name, std::ios::out);
  }
  void InitAnalysis() {
    //CellList::SetMinCellLength(10);
    n_bins_1d_ = params_->polar_order_n_bins;
    // Ensure n_bins_1d_ is even, to avoid headaches
    if (n_bins_1d_ % 2 != 0) {
      n_bins_1d_++;
    }
    n_bins_ = SQR(n_bins_1d_);
    histogram_ = new int[n_bins_];
    std::fill(histogram_, histogram_ + n_bins_, 0.0);
    contact_cut_ = params_->polar_order_contact_cutoff;
    contact_bin_width_ = contact_cut_ / n_bins_1d_;
    polar_bin_width_ = 2.0 / n_bins_1d_;
    output_ << "2D histogram: x-axis: contact_number, y-axis: local_polar_order\n";
    output_avg_ << "time avg_polar_order avg_contact_number\n";
    RequireInteractionAnalysis();
  }
  void RunAnalysis() {
    std::vector<double> po;
    std::vector<double> cn;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      //it->CalcPolarOrder();
      it->GetPolarOrders(po);
      it->GetContactNumbers(cn);
    }
    if (po.size() != cn.size()) {
      Logger::Error(
          "Number of polar order parameters and contact numbers not equal");
    }
    double po_avg = 0;
    double cn_avg = 0;
    for (int i = 0; i < po.size(); ++i) {
      po_avg += po[i];
      cn_avg += cn[i];
      if (cn[i] > contact_cut_) {
        continue;
      }
      int x = (int)(floor(cn[i] / contact_bin_width_));
      int y = (int)(floor((po[i] + 1) / polar_bin_width_));
      if (y == n_bins_1d_)
        y = n_bins_1d_ - 1;
      if (x == n_bins_1d_)
        x = n_bins_1d_ - 1;
      if (y == -1)
        y = 0;
      if (x == -1)
        x = 0;
      if (y < 0 || x < 0 || y > n_bins_1d_ - 1 || x > n_bins_1d_ - 1) {
        Logger::Warning("Contact number: %2.5f, Polar order: %2.5f", cn[i],
                        po[i]);
        Logger::Error("Encountered values outside of expected range in "
                      "PolarOrderAnalysis");
      }
      histogram_[n_bins_1d_ * y + x]++;
    }
    po_avg /= po.size();
    cn_avg /= cn.size();
    output_avg_ << time_ << " " << po_avg << " " << cn_avg << "\n";
  }

  void EndAnalysis() {
    /* In order to avoid overcounting cases where there were no local
     * interactors to count for local polar order, I am going to smooth the bin
     * representing (0,0) in the histogram, by averaging vertically along the
     * y-axis */
    int avg_bin = (int)floor(
        0.5 * (histogram_[(n_bins_1d_ / 2 - 1) * n_bins_1d_] +
               histogram_[(n_bins_1d_ / 2 + 1) * n_bins_1d_]));
    histogram_[(n_bins_1d_ / 2) * n_bins_1d_] = avg_bin;
    for (int i = 0; i < n_bins_1d_; ++i) {
      for (int j = 0; j < n_bins_1d_; ++j) {
        output_
            << histogram_[(n_bins_1d_ - 1 - i) * n_bins_1d_ + j]
            << " ";
      }
      output_ << "\n";
    }
    if (output_avg_.is_open()) {
      output_avg_.close();
    }
  }
};

#endif // _SIMCORE_FILAMENT_POLAR_ORDER_ANALYSIS_H_
