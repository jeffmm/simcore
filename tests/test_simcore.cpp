#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <simcore.hpp>

class Tester {
  public:
    static bool TestMngrInit(SimulationManager & mngr, run_options run_opts) {
      try {
        mngr.InitManager(run_opts);
      }
      catch (int e) {
        printf("Exception occurred with No. %d", e);
        return false;
      }
      return true;
    }
    void TestMngrParams(SimulationManager & mngr) {
      REQUIRE(mngr.run_name_.compare("test") == 0);
      REQUIRE(mngr.n_runs_ == 13);
      REQUIRE(mngr.n_random_ == 7);
      // We have two parameters for n_dim
      REQUIRE(mngr.n_var_ == 7*2);
    }
    static bool TestSim(Simulation & sim, system_parameters params) {
      sim.params_ = params;
      sim.run_name_ = params.run_name;
      try {
        sim.InitSimulation();
      }
      catch (int e) {
        printf("Exception occurred with No. %d", e);
        return false;
      }
      return true;
    }
};

TEST_CASE("Simulation manager") {
  run_options run_opts;
  run_opts.param_file = "../../tests/test_params.yaml";
  run_opts.default_param_file = "../../config/default_config.yaml";
  SimulationManager mngr;
  Tester test;
  SECTION("Simulation manager initialization") {
  REQUIRE(test.TestMngrInit(mngr, run_opts));
  test.TestMngrParams(mngr);
  }
}

TEST_CASE("Simulation") {
  system_parameters params;
  params.run_name = "test";
  Simulation sim;
  Tester test;
  SECTION("Simulation initialization") {
    REQUIRE(test.TestSim(sim, params));
  }
}


