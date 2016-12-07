#ifndef _SIMCORE_TEST_MANAGER_H_
#define _SIMCORE_TEST_MANAGER_H_

#include "auxiliary.h"
#include "helpers.h"
#include "test_module_base.h"

class TestManager {
  public:
    TestManager() {}
    ~TestManager() {}
   
    void InitManager(const std::string &filename);
    void RegisterTestModules();
    void RunTestModules();

  private:
    long seed_;
    std::string filename_;

    rfh::factory test_module_factory_;
    std::vector<TestModuleBase*> test_modules_;
};

#endif /* SIMCORE_TEST_MANAGER_H_ */
