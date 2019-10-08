#ifndef _SIMCORE_EXPONENTIAL_DIST_H_
#define _SIMCORE_EXPONENTIAL_DIST_H_

#include "auxiliary.hpp"

/* The exponential distribution has the functional form
     PDF(x) = \lambda e^{-\lambda x},
   where \lambda is a "polydispersity factor." It has a cumulative distribution
   function
     CDF(x) = 1 - e^{-\lambda x}
   This class is meant to sample randomly from the exponential distribution by
   inverting the CDF so that it returns
     CDF^{-1}( CDF(u) ) == u
   where u is a uniformly distributed random number between 0 and 1. Since
   the CDF is analytically invertible, this simply returns:
     CDF^{-1}(u) = -1/\lambda ln(1 - u) = x

*/
class ExponentialDist {
private:
  double inv_lambda_ = 50; // polydispersity parameter, mean of distribution
  double min_roll_ = 0; // smallest allowed roll for lowest allowed value of x
  double max_roll_ = 1; // largest allowed roll for largest allowed value of x
  double CDF(double x) { return 1 - exp(-1.0 / inv_lambda_ * x); }
  double InvCDF(double x) { return -inv_lambda_ * log(1 - x); }

public:
  /* Initialize the exponential distribution generator. The only necessary
     parameter is the polydispersity parameter, which is the mean of the
     distribution (1/lambda) */
  void Init(double mean, double x_min = 0, double x_max = 1e6) {
    inv_lambda_ = mean;
    if (mean < 1 || mean > 10000) {
      Logger::Error(
          "ExponentialDist given off-base polydispersity parameter (%2.2f). "
          "If you /really/ want this parameter, you'll have to remove this "
          "error",
          mean);
    }
    if (x_min < 0) {
      Logger::Error("ExponentialDist expects a non-negative minimum x value");
    } else if (x_min > x_max) {
      Logger::Error("Value for x_min in ExponentialDist is larger than x_max");
    } else if (x_min > 0) {
      min_roll_ = CDF(x_min);
      max_roll_ = CDF(x_max);
      Logger::Trace("Renormalizing rolls in ExponentialDist generator to have "
                    "a minimum of %2.2f and maximum of %2.2f to avoid lengths "
                    "shorter than minimum length of %2.2f and maximum length "
                    "of %2.2f",
                    min_roll_, max_roll_, x_min, x_max);
    }
  }
  /* Given a uniform random number in the range of 0 and 1 (roll), return its
     corresponding exponentially-distributed random value */
  double Rand(double roll) {
    /* If we reject values less than some number (because of a minimum length
       requirement) then renormalize the roll to be within the acceptable range.
       This is done for both the lower and upper bound on length */
    if (min_roll_ > 0) {
      roll = (max_roll_ - min_roll_) * roll + min_roll_;
    }
    double result = InvCDF(roll);
    Logger::Trace("ExponentialDist generator received a roll of %2.2f and "
                  "returned a length of %2.2f",
                  roll, result);
    return result;
  }
};

#endif
