#ifndef _SIMCORE_TEST_XLINK_KMC_H_
#define _SIMCORE_TEST_XLINK_KMC_H_

#include "xlink_kmc.h"
#include "test_base.h"

class TestXlinkKMC : public TestBase, public XlinkKMC {
  public:
    TestXlinkKMC() : TestBase(), XlinkKMC() {}

    virtual void InitTests();
    virtual void UnitTests();
    virtual void IntegrationTests();

  protected:

};

#endif /* _SIMCORE_TEST_XLINK_KMC_H_ */
