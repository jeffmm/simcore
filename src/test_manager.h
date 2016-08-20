#ifndef _SIMCORE_TEST_MANAGER_H_
#define _SIMCORE_TEST_MANAGER_H_

#include "auxiliary.h"
#include "helpers.h"
#include "test_module_base.h"

class TestManager {
  public:
    TestManager() {}
    ~TestManager() {
      std::cout << "TestManager destructor\n";
      for (auto i = 0; i < test_modules_.size(); ++i) {
        delete(test_modules_[i]);
      }
      test_modules_.clear();
    }
   
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
