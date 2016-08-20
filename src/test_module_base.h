#ifndef _SIMCORE_TEST_MODULE_BASE_H_
#define _SIMCORE_TEST_MODULE_BASE_H_

#include "auxiliary.h"

#include <yaml-cpp/yaml.h>

class TestModuleBase {
  protected:

    std::string name_;
    YAML::Node node_;

  public:
    TestModuleBase() {}
    virtual ~TestModuleBase() {}

    virtual void InitTestModule(const std::string& filename) = 0;
    virtual void RunTests() = 0;
    virtual void UnitTests() = 0;
    virtual void IntegrationTests() = 0;
};

#endif /* _SIMCORE_TEST_MODULE_BASE_H_ */
