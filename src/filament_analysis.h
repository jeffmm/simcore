#ifndef _SIMCORE_FILAMENT_ANALYSIS_H_
#define _SIMCORE_FILAMENT_ANALYSIS_H_

#include "posit_reader.h"

class FilamentAnalysis {

  private:
    PositReader preader_;
    std::string posit_file_name_;
    int n_dim_;
    int n_objs_;
    int n_posit_;
    int n_steps_;
    int n_time_;
    int n_validate_;
    int n_bins_;
    int time_;
    double delta_;
    double plength_;
    double * diameters_ = nullptr;
    double * lengths_;
    double * child_lengths_;
    double * sqr_end_to_end_;
    double * cos_thetas_;
    double * misc_;
    int * theta_distribution_;
    double ms_end_to_end_;
    double ms_end_to_end_err_;
    double gamma_par_;
    double gamma_perp_;
    double gamma_rot_;

    void CleanUp() {
      preader_.CloseFile();
      delete[] diameters_;
      delete[] lengths_;
      delete[] child_lengths_;
      delete[] sqr_end_to_end_;
      delete[] cos_thetas_;
      delete[] misc_;
      delete[] theta_distribution_;
      diameters_ = nullptr;
    }
    // Returns index of contiguous array from 2d coords
    void FilamentInit() {
      n_steps_ = preader_.NSteps();
      n_posit_ = preader_.NPosit();
      n_time_ = n_steps_/n_posit_;
      // Make use of integer division here
      preader_.GetNObjs(&n_objs_);
      if (diameters_ == nullptr) {
        diameters_ = new double[n_objs_];
        lengths_ = new double[n_objs_];
        misc_ = new double[n_objs_*3];
        child_lengths_ = new double[n_objs_];
        sqr_end_to_end_ = new double[n_objs_];
        cos_thetas_ = new double[n_objs_];
        theta_distribution_ = new int[n_bins_];
        std::fill(theta_distribution_, theta_distribution_+n_bins_, 0);
      }
    }

    // Put the info from the misc array into their relative arrays
    // The order is child_length, end_to_end, cos_theta
    void OrganizeMisc() {
      for (int i=0; i<n_objs_; ++i) {
        child_lengths_[i] = misc_[3*i];
        sqr_end_to_end_[i] = misc_[3*i+1];
        cos_thetas_[i] = misc_[3*i+2];
      }
    }

    void BinTheta() {
      for (int i=0; i<n_objs_; ++i) {
        int bin_number = (int) floor( (1 + cos_thetas_[i]) * (n_bins_/2) );
        // Check boundaries
        if (bin_number == n_bins_)
          bin_number = n_bins_-1;
        else if (bin_number == -1)
          bin_number = 0;
        // Check for nonsensical values
        else if (bin_number > n_bins_ || bin_number < 0) error_exit("Something went wrong in BinTheta!\n");
        theta_distribution_[bin_number]++;
      }
    }
   // Calculate mean square distance <(r(t)-r(o))^2> and standard error
    void CalcMeanSquareEndToEnd() {
      for (int i=0; i<n_objs_; ++i) {
        ms_end_to_end_ += sqr_end_to_end_[i];
        ms_end_to_end_err_ += SQR(sqr_end_to_end_[i]);
      }
    }
    void WriteData() {
      std::ostringstream file_name;
      file_name << posit_file_name_.substr(0,posit_file_name_.find_last_of(".")) << ".filament";
      std::ofstream fil_file(file_name.str().c_str(), std::ios_base::out);
      fil_file << "# n_dim delta n_steps n_posit n_objs \n";
      fil_file << n_dim_ << " " << delta_ << " " << n_steps_ << " " << n_posit_ << " " << n_objs_ << " " << "\n";
      fil_file << "# mean_square_end_to_end ms_end_to_end_err length child_length persistence_length\n";
      fil_file << ms_end_to_end_ << " " << ms_end_to_end_err_ << " " << lengths_[0] << " " << child_lengths_[0] << " "<< plength_ << "\n";
      fil_file << "# theta_distribution_histogram\n";
      for (int i=0; i<n_bins_; ++i)
        fil_file << theta_distribution_[i] << "\n";
      fil_file.close();
    }
    void FinalizeData() {
      ms_end_to_end_/=time_;
      ms_end_to_end_/=n_objs_;
      ms_end_to_end_err_/=time_;
      ms_end_to_end_err_/=n_objs_;
      ms_end_to_end_err_ = sqrt((ms_end_to_end_err_ - SQR(ms_end_to_end_))/n_objs_);
    }

  public:
    void AnalyzeFilament(system_parameters *params, std::string posit_file_name) {
      delta_ = params->delta;
      n_dim_ = params->n_dim;
      n_bins_ = params->n_bins;
      posit_file_name_ = posit_file_name;
      plength_ = params->persistence_length;
      preader_.LoadFile(posit_file_name);
      SID posit_sid = preader_.GetSID();
      if (posit_sid != SID::filament) {
        std::cout << "      Aborting filament analysis on " << posit_file_name << "\n";
        preader_.CloseFile();
        return;
      }
      FilamentInit();
      time_=0;
      int nobj;
      while (preader_.ReadNext()) {
        preader_.GetNObjs(&nobj);
        if (nobj != n_objs_) {
          printf("n_objs_new: %d, n_objs_prev: %d\n",nobj, n_objs_);
          error_exit("ERROR: Number of objects changed in diffusion analysis!\n");
        }
        else if (time_ > n_time_)
          error_exit("ERROR: Didn't hit EOF when expected in diffusion analysis.\n");
        else {
          preader_.GetDiameter(diameters_);
          preader_.GetLength(lengths_);
          preader_.GetMisc(misc_);
          OrganizeMisc();
          CalcMeanSquareEndToEnd();
          BinTheta();
          time_++;
        }
      }
      FinalizeData();
      WriteData();
      CleanUp();
    }
};

#endif // _SIMCORE_FILAMENT_ANALYSIS_H_

