#include "test_xlink_kmc.h"

#include "test_helpers.h"

#include <iomanip>
#include <yaml-cpp/yaml.h>

void TestXlinkKMC::InitTestModule(const std::string& filename) {
  name_ = "TestXlinkKMC";
  std::cout << name_ << " Init Tests\n";

  // Add the tests outlined in the yaml file?
  node_ = YAML::LoadFile(filename);

  // Test out functional
  if (node_[name_]["CalcCutoff"]) {
    std::function<bool(void)> f_calc_cutoff = std::bind(&TestXlinkKMC::UnitTestCalcCutoff, this);
    unit_tests_.push_back(f_calc_cutoff);
  }
}

void TestXlinkKMC::RunTests() {
  std::cout << name_ << " Run Tests\n";
  UnitTests();
  IntegrationTests();
}

void TestXlinkKMC::UnitTests() {
  std::cout << name_ << " Unit Tests ->\n";
  bool test_success = true;
  for (auto fit = unit_tests_.begin(); fit != unit_tests_.end(); ++fit) {
    test_success |= (*fit)();
  }
  std::cout << name_;
  std::cout << " Unit Tests: " << (test_success ? "passed" : "failed") << std::endl;
}

void TestXlinkKMC::IntegrationTests() {
  std::cout << name_;
  std::cout << " Integration Tests\n";
}

bool TestXlinkKMC::UnitTestCalcCutoff() {
  std::cout << "  Unit Testing XlinkKMC CalcCutoff\n";
  bool success = true;

  // Grab the YAML file information
  eps_eff_1_2_[0] = eps_eff_1_2_[1] = node_[name_]["CalcCutoff"]["concentration_1_2"].as<double>();
  k_stretch_  = node_[name_]["CalcCutoff"]["spring_constant"].as<double>();
  max_length_ = node_[name_]["CalcCutoff"]["max_length"].as<double>();

  std::vector<double> equilibrium_lengths;
  std::vector<double> barrier_weights;
  std::vector<double> results;

  for (int i = 0; i < node_[name_]["CalcCutoff"]["results"].size(); ++i) {
    barrier_weights.push_back(node_[name_]["CalcCutoff"]["results"][i][0].as<double>());
    equilibrium_lengths.push_back(node_[name_]["CalcCutoff"]["results"][i][1].as<double>());
    results.push_back(node_[name_]["CalcCutoff"]["results"][i][2].as<double>());
  }

  for (int i = 0; i < results.size(); ++i) {
    barrier_weight_ = barrier_weights[i];
    r_equil_ = equilibrium_lengths[i];

    CalcCutoff();

    double res = rcutoff_1_2_;

    if (uth::almost_equal<double>(res, results[i], 1e-8)) {
      std::cout << "    Test[" << i << "] passed\n";
    } else {
      std::cout << "    Test[" << i << "] failed\n";
      success = false;
    }
  }

  return success;
}
