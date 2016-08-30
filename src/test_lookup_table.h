#ifndef _SIMCORE_TEST_LOOKUP_TABLE_H_
#define _SIMCORE_TEST_LOOKUP_TABLE_H_

#include "lookup_table.h"
#include "space.h"
#include "test_module_base.h"
#include "xlink_helpers.h"

#include <functional>

class TestLookupTable : public TestModuleBase, public LookupTable {
  public:
    TestLookupTable() : TestModuleBase(), LookupTable() {}
    virtual ~TestLookupTable() {}

    virtual void InitTestModule(const std::string& filename);
    virtual void RunTests();
    virtual void UnitTests();
    virtual void IntegrationTests();

  protected:

    bool UnitTestBobFacsimile(int test_num);

    std::vector<std::function<bool(int)>> unit_tests_;
    std::vector<std::string> unit_tests_names_;
    std::vector<std::vector<bool>> unit_tests_results_;
};

#endif
