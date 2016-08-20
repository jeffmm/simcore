#ifndef _SIMCORE_TEST_XLINK_KMC_H_
#define _SIMCORE_TEST_XLINK_KMC_H_

#include "xlink_kmc.h"
#include "test_module_base.h"

#include <functional>

class TestXlinkKMC : public TestModuleBase, public XlinkKMC {
  public:
    TestXlinkKMC() : TestModuleBase(), XlinkKMC() {}
    virtual ~TestXlinkKMC() {}

    virtual void InitTestModule(const std::string& filename);
    virtual void RunTests();
    virtual void UnitTests();
    virtual void IntegrationTests();

  protected:

    bool UnitTestCalcCutoff();

    std::vector<std::function<bool(void)>> unit_tests_;
};

#endif /* _SIMCORE_TEST_XLINK_KMC_H_ */
