#include "test_xlink_kmc.h"

#include "br_rod.h"

#include <iomanip>
#include <yaml-cpp/yaml.h>

void TestXlinkKMC::InitTests(const std::string& filename) {
  std::cout << "TestXlinkKMC Init Tests\n";
  // Spoof the params file to what we want
}

void TestXlinkKMC::UnitTests() {
  std::cout << "TestXlinkKMC Unit Tests ->\n";
  bool test_success = true;
  test_success = UnitTestCalcCutoff();
  std::cout << "TestXlinkKMC Unit Tests: " << (test_success ? "passed" : "failed") << std::endl;
}

void TestXlinkKMC::IntegrationTests() {
  std::cout << "TestXlinkKMC Integration Tests\n";
}

bool TestXlinkKMC::UnitTestCalcCutoff() {
  std::cout << "  Unit Testing XlinkKMC CalcCutoff\n";
  bool success = true;

  // Prep everything we need for calc cutoff
  eps_eff_1_2_[0] = eps_eff_1_2_[1] = 10.0;
  k_stretch_ = 45.29;
  max_length_ = 110;

  // barrier weight, requil
  std::map<double, std::map<double, double>> results_map;
  results_map[0.0][0.0] = 0.758961;
  results_map[0.0][3.12] = 3.87896;
  results_map[0.0001][0.0] = 0.758999;
  results_map[0.0001][3.12] = 3.879;
  results_map[0.25][0.0] = 0.876372;
  results_map[0.25][3.12] = 3.99637;

  int itest = 0;
  for (auto bit = results_map.begin(); bit != results_map.end(); ++bit) {
    for (auto rit = bit->second.begin(); rit != bit->second.end(); ++rit) {
      barrier_weight_ = bit->first;
      r_equil_ = rit->first;

      CalcCutoff();

      double result = rcutoff_1_2_;

      if (uth::almost_equal<double>(result, rit->second, 1E-3)) {
        std::cout << "    Test[" << itest << "] passed\n";
      } else {
        std::cout << "    Test[" << itest << "] failed\n";
        success = false;
      }
      itest++;
    }
  }

  return success;
}
