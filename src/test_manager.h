#ifndef _SIMCORE_TEST_MANAGER_H_
#define _SIMCORE_TEST_MANAGER_H_

#include "auxiliary.h"
#include "helpers.h"

class TestManager {
  public:
    TestManager() {}
   
    void InitManager(const std::string &filename);
    void RegisterTests();
    void RunTests();

  private:
    long seed_;

    rfh::factory test_factory_;
};

#endif /* SIMCORE_TEST_MANAGER_H_ */
