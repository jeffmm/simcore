#ifndef _SIMCORE_TEST_HELPERS_H_
#define _SIMCORE_TEST_HELPERS_H_

// namespace unit test helpers
namespace uth {
  inline bool compare_doubles(double x, double y) {
    return (x >= y - std::numeric_limits<double>::epsilon() &&
            x <= y + std::numeric_limits<double>::epsilon()) ? true : false;
  }

  template<typename T>
  inline bool almost_equal(T x, T y, T eps) {
    return (x >= y - eps && x <= y + eps) ? true : false;
  }
}


#endif
