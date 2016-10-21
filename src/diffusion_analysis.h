#ifndef _SIMCORE_DIFFUSION_ANALYSIS_H_
#define _SIMCORE_DIFFUSION_ANALYSIS_H_

#include "posit_reader.h"

class DiffusionAnalysis {

  private:
    PositReader preader_;
    std::string posit_file_name_;
    int n_dim_;
    int n_objs_;
    int n_posit_;
    int n_steps_;
    int n_time_;
    int n_validate_;
    int n_interval_; // Number of time avgs depending on n_interval and n_steps
    int time_avg_interval_; // The number of iterations over which to time average
    int time_;
    double delta_;
    double * positions_0_ = nullptr;
    double * orientations_0_;
    double * positions_;
    double * scaled_pos_;
    double * orientations_;
    double * diameters_;
    double * lengths_;
    double * msd_; // mean square displacement
    double * vcf_; // vector correlation function
    double * msd_err_;
    double * vcf_err_;
    double gamma_par_;
    double gamma_perp_;
    double gamma_rot_;

    void SetInitPos() {
      preader_.GetNObjs(&n_objs_);
      if (positions_0_ == nullptr) {
        positions_0_ = new double[n_objs_*3];
        orientations_0_ = new double[n_objs_*3];
        positions_ = new double[n_objs_*3];
        scaled_pos_ = new double[n_objs_*3];
        orientations_ = new double[n_objs_*3];
        diameters_ = new double[n_objs_];
        lengths_ = new double[n_objs_];
      }
      preader_.GetPosit(positions_, scaled_pos_, orientations_, diameters_, lengths_);
      std::copy(positions_, positions_+3*n_objs_, positions_0_);
      std::copy(orientations_, orientations_+3*n_objs_, orientations_0_);
    }
    void CleanUp() {
      preader_.CloseFile();
      delete[] positions_0_;
      delete[] orientations_0_;
      delete[] positions_;
      delete[] orientations_;
      delete[] scaled_pos_;
      delete[] diameters_;
      delete[] lengths_;
      delete[] vcf_;
      delete[] vcf_err_;
      delete[] msd_;
      delete[] msd_err_;
    }
    // Returns index of contiguous array from 2d coords
    int ix(int x, int y) {
      return 3*x+y;
    }
    void DiffusionInit() {
      n_steps_ = preader_.NSteps();
      n_posit_ = preader_.NPosit();
      n_time_ = n_steps_/n_posit_;
      // Make use of integer division here
      n_interval_ = n_time_/time_avg_interval_;
      n_interval_ = (n_interval_ < 1 ? 1 : n_interval_);
      time_avg_interval_ = n_time_/n_interval_;
      SetInitPos();
      AllocateAnalysis();
    }
    void AllocateAnalysis() {
      vcf_ = new double[time_avg_interval_];
      msd_ = new double[time_avg_interval_];
      vcf_err_ = new double[time_avg_interval_];
      msd_err_ = new double[time_avg_interval_];
      std::fill(msd_,msd_+time_avg_interval_,0.0);
      std::fill(vcf_,vcf_+time_avg_interval_,0.0);
      std::fill(msd_err_,msd_err_+time_avg_interval_,0.0);
      std::fill(vcf_err_,vcf_err_+time_avg_interval_,0.0);
    }
    // Calculate vector correlation function <u(t).u(o)> and standard error
    void CalcVCF() {
      double avg_udotu0 = 0.0;
      double avg_udotu0_sqr = 0.0;
      for (int i=0; i<n_objs_; ++i) {
        double udotu0 = 0.0;
        for (int j=0; j<n_dim_; ++j) {
          udotu0 += orientations_[ix(i,j)] * orientations_0_[ix(i,j)];
        }
        avg_udotu0 += udotu0;
        avg_udotu0_sqr += udotu0*udotu0;
      }
      avg_udotu0/=n_objs_;
      avg_udotu0_sqr/=n_objs_;
      double stdev2 = avg_udotu0_sqr - SQR(avg_udotu0);
      vcf_[time_] += avg_udotu0/stdev2;
      vcf_err_[time_] += 1.0/stdev2;
    }
    // Calculate mean square distance <(r(t)-r(o))^2> and standard error
    void CalcMSD() {
      double avg_sqr_dist = 0.0;
      double avg_sqr_dist_sqr = 0.0;
      for (int i=0; i<n_objs_; ++i) {
        for (int j=0; j<n_dim_; ++j) {
          double r_diff = positions_[ix(i,j)] - positions_0_[ix(i,j)];
          r_diff = SQR(r_diff);
          avg_sqr_dist += r_diff;
          avg_sqr_dist_sqr += SQR(r_diff);
        }
      }
      avg_sqr_dist/=n_objs_;
      avg_sqr_dist_sqr/=n_objs_;
      double stdev2 = avg_sqr_dist_sqr - SQR(avg_sqr_dist);
      msd_[time_] += avg_sqr_dist/stdev2;
      msd_err_[time_] += 1.0/stdev2;
    } 
    void WriteData() {
      std::ostringstream file_name;
      file_name << posit_file_name_ << ".diffusion";
      std::ofstream diff_file(file_name.str().c_str(), std::ios_base::out);
      diff_file << "# n_dim delta n_steps n_posit n_objs n_interval\n";
      diff_file << n_dim_ << " " << delta_ << " " << n_steps_ << " " << n_posit_ << " " << n_objs_ << " " << n_interval_ << "\n";
      diff_file << "# time msd msd_err vcf vcf_err\n";
      diff_file << "0.0 0.0 0.0 1.0 0.0\n";
      for (int t=0; t<time_avg_interval_; ++t)
        diff_file << (t+1)*delta_*n_posit_ << " " << msd_[t] << " " << msd_err_[t]
          << " " << vcf_[t] << " " << vcf_err_[t] << "\n";
      diff_file.close();
    }
    void FinalizeData() {
      // vcf_err_ is 1/std^2 = sum_i{1/std_i^2}
      // vcf_ is sum_i{a_i / std^2}
      // need vcf_ to be a = sum_i{a_i / std^2} / sum_i {1.0/std^2}
      // need vcf_err_ to be SEM = std/sqrt(n_objs)
      for (int t=0; t<time_avg_interval_; ++t) {
        msd_[t] = msd_[t]/msd_err_[t];
        vcf_[t] = vcf_[t]/vcf_err_[t];
        msd_err_[t] = 1.0/sqrt(n_objs_*msd_err_[t]);
        vcf_err_[t] = 1.0/sqrt(n_objs_*vcf_err_[t]);
      }
    }

  public:
    void CalculateDiffusion(system_parameters *params, std::string posit_file_name) {
      delta_ = params->delta;
      n_dim_ = params->n_dim;
      time_avg_interval_ = params->diffusion_interval;
      posit_file_name_ = posit_file_name;
      preader_.LoadFile(posit_file_name);
      DiffusionInit();
      time_=0;
      int nobj;
      while (preader_.GetNext(&nobj, positions_,
             scaled_pos_, orientations_, diameters_, lengths_)) {
        if (nobj != n_objs_) {
          printf("n_objs_new: %d, n_objs_prev: %d\n",nobj, n_objs_);
          error_exit("ERROR: Number of objects changed in diffusion analysis!\n");
        }
        else if (time_ > n_time_ || time_ > time_avg_interval_)
          error_exit("ERROR: Didn't hit EOF when expected in diffusion analysis.\n");
        if (time_ == time_avg_interval_) {
          if (n_interval_ == 1)
            error_exit("Something went wrong in diffusion analysis\n");
          SetInitPos();
          time_=0;
        }
        else {
          CalcVCF();
          CalcMSD();
          time_++;
        }
      }
      FinalizeData();
      WriteData();
      CleanUp();
    }
};

#endif // _SIMCORE_DIFFUSION_ANALYSIS_H_

