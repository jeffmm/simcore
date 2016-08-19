#ifndef _SIMCORE_TEST_MANAGER_H_
#define _SIMCORE_TEST_MANAGER_H_

#include "auxiliary.h"
#include "helpers.h"
#include "test_base.h"

class TestManager {
  public:
    TestManager() {}
   
    void InitManager(const std::string &filename);
    void RegisterTests();
    void RunTests();

  private:
    long seed_;
    std::string filename_;

    rfh::factory test_factory_;
    std::vector<TestBase*> tests_;
};

#endif /* SIMCORE_TEST_MANAGER_H_ */
