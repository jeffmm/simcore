 
#include <iostream>
#include <stdlib.h>
#include "simulation_manager.h"

/*************************
   ::SimCORE Main::
   Parse commandline flags and start simulation manager
**************************/
template <class TRod>
TRod MockRod(int id) {
  TRod rod;
  rod.gid = id;
  rod.length = 400;
  rod.rank = 0;
  rod.radius = .5;
  for (int i = 0; i < 3; ++i) {
    rod.pos[i] = 0;
    rod.direction[i] = 0;
  }
  rod.direction[0] = 1;
  return rod;
}

void test_kmcx() {
  // Time step
  double dt = .0001;
  double KBT = 1.;
  // Protein data
  ExampleXlink xlink;
  xlink.setMockXlink();
  // Sylinder data
  ExampleRod rod = MockRod<ExampleRod>(0);
  KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
  double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
  printf("%d %d %d\n", (int) (kmc.getMu(0) == 0), (int) (kmc.getDistMin(0) == 0), (int) (ABS(prob - 0.0001909859317102744) <= 1e-8));
}

int main(int argc, char *argv[]) {

  // Parse input flags, see parse_flags.h for documentation
  run_options run_opts = parse_opts(argc, argv);

  // Initialize manager assets
  SimulationManager sim;
  sim.InitManager(run_opts);

  // Main control function
  sim.RunManager();

  test_kmcx();
  return 0;
}

