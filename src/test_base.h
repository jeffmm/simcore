#ifndef _SIMCORE_TEST_BASE_H_
#define _SIMCORE_TEST_BASE_H_

// Base test class to inherit from for everything, yay
class TestBase {
  protected:

  public:
    TestBase() {}

    virtual void InitTests(const std::string& filename) = 0;
    virtual void UnitTests() = 0;
    virtual void IntegrationTests() = 0;
};

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

#endif /* _SIMCORE_TEST_BASE_H_ */
