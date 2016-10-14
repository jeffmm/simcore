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
    int time_avg_interval_; // The number of iterations over which to time average
    int time_;
    double delta_;
    double * positions_0_;
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
      n_dim_ = 2;
      if (!preader_.GetNObjs(&n_objs_)) {
        printf("Something went wrong in diffusion init.\n");
        exit(1);
      }
      positions_0_ = new double[n_objs_*3];
      orientations_0_ = new double[n_objs_*3];
      positions_ = new double[n_objs_*3];
      scaled_pos_ = new double[n_objs_*3];
      orientations_ = new double[n_objs_*3];
      diameters_ = new double[n_objs_];
      lengths_ = new double[n_objs_];
      if (!preader_.GetPosit(positions_, scaled_pos_, orientations_, diameters_, lengths_)) {
        printf("Something went wrong in diffusion init.\n");
        exit(1);
      }
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
      SetInitPos();
      AllocateAnalysis();
    }
    void AllocateAnalysis() {
      vcf_ = new double[n_time_];
      msd_ = new double[n_time_];
      vcf_err_ = new double[n_time_];
      msd_err_ = new double[n_time_];
    }
    void CalcVCF() {
      double inv_sqrt_nobj = 1.0/sqrt(n_objs_);
      double avg_udotu0 = 0.0;
      double avg_udotu0_sqr = 0.0;
      for (int i=0; i<n_objs_; ++i) {
        double udotu0 = 0.0;
        for (int j=0; j<n_dim_; ++j) {
          udotu0 += orientations_[ix(i,j)] * orientations_0_[ix(i,j)];
          //printf("%2.8f ",orientations_0_[ix(i,j)]);
        }
        //printf("\n");
        avg_udotu0 += udotu0;
        avg_udotu0_sqr += udotu0*udotu0;
      }
      //exit(1);
      avg_udotu0/=n_objs_;
      avg_udotu0_sqr/=n_objs_;
      vcf_[time_] = avg_udotu0;
      vcf_err_[time_] = inv_sqrt_nobj * sqrt(avg_udotu0_sqr - SQR(avg_udotu0));
    }
    void CalcMSD() {} // TODO
    void WriteData() {
      std::ostringstream file_name;
      file_name << posit_file_name_ << ".diffusion";
      std::ofstream diff_file(file_name.str().c_str(), std::ios_base::out);
      diff_file << n_objs_ << " " << n_steps_ << " " << n_posit_ << "\n";
      for (int t=0; t<n_time_; ++t)
        diff_file << vcf_[t] << " " << vcf_err_[t] << "\n";
      diff_file.close();
    }

  public:
    void CalculateDiffusion(std::string posit_file_name) {
      posit_file_name_ = posit_file_name;
      preader_.LoadFile(posit_file_name);
      DiffusionInit();
      time_=0;
      int nobj;
      while (preader_.GetNext(&nobj, positions_,
             scaled_pos_, orientations_, diameters_, lengths_)) {
        if (nobj != n_objs_) {
          printf("%d %d\n",nobj, n_objs_);
          printf("ERROR: Number of objects changed in diffusion analysis!\n");
          exit(1);
        }
        else if (time_ > n_time_) {
          printf("ERROR: Didn't hit EOF when expected in diffusion analysis!\n");
          exit(1);
        }
        CalcVCF();
        CalcMSD();
        time_++;
      }
      WriteData();
      CleanUp();
    }
};

#endif // _SIMCORE_DIFFUSION_ANALYSIS_H_

