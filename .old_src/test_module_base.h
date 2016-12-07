#ifndef _SIMCORE_TEST_MODULE_BASE_H_
#define _SIMCORE_TEST_MODULE_BASE_H_

#include "auxiliary.h"

#include <functional>

#include <yaml-cpp/yaml.h>

class TestModuleBase {
  protected:

    std::string name_;
    std::string filename_;
    YAML::Node node_;

    std::vector<std::function<bool(int)>> tests_;
    std::vector<std::string> test_names_;
    std::vector<std::vector<bool>> test_results_;

  public:
    TestModuleBase() {}
    virtual ~TestModuleBase() {}

    virtual void InitTestModule(const std::string& filename) = 0;
    virtual void RunTests() = 0;
};

#endif /* _SIMCORE_TEST_MODULE_BASE_H_ */
