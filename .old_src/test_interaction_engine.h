#ifndef _SIMCORE_TEST_INTERACTION_ENGINE_H_
#define _SIMCORE_TEST_INTERACTION_ENGINE_H_

#include "interaction_engine.h"
#include "space.h"
#include "test_module_base.h"

#include <functional>

class BrRod;
class SpindlePoleBody;

class TestInteractionEngine : public TestModuleBase, public InteractionEngine {
  public:
      TestInteractionEngine() : TestModuleBase(), InteractionEngine() {}
      virtual ~TestInteractionEngine() {}

      virtual void InitTestModule(const std::string& filename);
      virtual void RunTests();

    protected:

      bool TestInit(int test_num);
      bool TestTetherParticlesMP(int test_num);

      system_parameters params_sub_;
      SpaceProperties space_sub_;
};

#endif
