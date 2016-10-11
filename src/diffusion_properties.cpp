#include "auxiliary.h"

class DiffusionProperties {

  private:
    int n_dim_;
    int n_objects_;
    int n_validate_;
    int time_avg_interval_; // The number of iterations over which to time average
    double delta_;
    double ** positions_0_;
    double ** orientations_0_;
    double ** positions_;
    double ** orientations_;
    double * msd_; // mean square displacement
    double * vcf_; // vector correlation function
    double gamma_par_;
    double gamma_perp_;
    double gamma_rot_;

  public:


};
