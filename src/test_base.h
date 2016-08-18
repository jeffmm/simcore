#ifndef _SIMCORE_TEST_BASE_H_
#define _SIMCORE_TEST_BASE_H_

// Base test class to inherit from for everything, yay

class TestBase {
  protected:

  public:
    TestBase() {
    }

    virtual void InitTests() = 0;
    virtual void UnitTests() = 0;
    virtual void IntegrationTests() = 0;
};

#endif /* _SIMCORE_TEST_BASE_H_ */
