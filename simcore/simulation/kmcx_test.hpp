#include <kmcx.hpp>

/* TESTS FOR KMCX-SIMCORE INTEGRATION */
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
  double dt = .0001;
  double KBT = 1.;
  ExampleXlink xlink;
  xlink.setMockXlink();
  ExampleRod rod = MockRod<ExampleRod>(0);
  KMC<ExampleRod> kmc(xlink.getPosPtr(), 1, xlink.getRcutUS());
  double prob = kmc.CalcProbUS(0, rod, xlink.getBindingFactorUS(0, dt));
  int test1 = (int)(kmc.getMu(0) == 0);
  int test2 = (int)(kmc.getDistMin(0) == 0);
  int test3 = (int)(ABS(prob - 0.0001909859317102744) <= 1e-8);

  if (test1 + test2 + test3 == 3) {
    printf("  All tests for KMCX-Simcore integration have passed!\n");
  } else {
    printf("  Tests for KMCX-Simcore integration have failed!\n");
  }
}
