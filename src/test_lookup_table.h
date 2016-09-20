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

  protected:

    bool UnitTestBobFacsimile(int test_num);
};

#endif
