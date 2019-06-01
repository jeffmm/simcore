#define CATCH_CONFIG_MAIN
#include <catch.hpp>


#include "../src/simcore/simulation_manager.h"
TEST_CASE("Simulation manager") {
  SECTION("Simulation manager initialization") {
    run_options run_opts;
    run_opts.param_file = "../test_params.yaml";
    SimulationManager mngr;
    mngr.InitManager(run_opts);
  }
}

