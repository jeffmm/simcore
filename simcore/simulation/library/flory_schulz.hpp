#include <math.h>

/* The Flory-Schulz distribution is the distribution used for polydisperse
   polymers.

   w_a(k) = a^2 k (1-a)^(k-1)

   CDF = 1 - (1-a)^k (1+ak)


*/
class FlorySchulz {
private:
  double a_ = 0.5;
  double k0_ = 50;
  double epsilon_ = 1e-4;
  double CDF(double k) { return 1 - pow(1 - a_, k) * (1 + a_ * k); }
  double Func(double k) { return a_ * a_ * k * pow(1 - a_, k - 1); }
  double NewtonRaphson(double roll) {
    double h;
    double k = k0_;
    do {
      h = (CDF(k) - roll) / Func(k);
      k -= h;
    } while (abs(h) >= epsilon_);
    return k;
  }

public:
  void Init(double a, double k0, double epsilon) {
    a_ = a;
    k0_ = k0;
    epsilon = epsilon_;
  }
  double Rand(double roll) { return NewtonRaphson(roll); }
};
