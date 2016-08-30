#include "test_lookup_table.h"

#include "test_helpers.h"

void TestLookupTable::InitTestModule(const std::string& filename) {
  name_ = "TestLookupTable";
  std::cout << filename << std::endl;
  std::cout << name_ << " Init Tests\n";

  // Add the tests outlined in the yaml file?
  node_ = YAML::LoadFile(filename);

  if (node_[name_]["BobFacsimile"]) {
    std::function<bool(int)> f_bobfacsimile = std::bind(&TestLookupTable::UnitTestBobFacsimile, this, std::placeholders::_1);
    unit_tests_.push_back(f_bobfacsimile);
    unit_tests_names_.push_back("BobFacsimile");
    unit_tests_results_.push_back(std::vector<bool>());
  }
}

void TestLookupTable::RunTests() {
  std::cout << name_ << " Run Tests\n";
  UnitTests();
  IntegrationTests();
}

void TestLookupTable::UnitTests() {
  std::cout << name_ << " Unit Tests ->\n";
  bool test_success = true;

  for (int i = 0; i < unit_tests_.size(); ++i) {
    bool this_success = unit_tests_[i](i);
    test_success &= this_success;

    std::cout << "----------------\n";
    std::cout << "--Unit Test " << unit_tests_names_[i] << ": ";
    std::cout << (this_success ? "PASSED" : "FAILED") << std::endl;
    for (int didpass = 0; didpass < unit_tests_results_[i].size(); ++didpass) {
      std::cout << "    Test[" << didpass << "] " << (unit_tests_results_[i][didpass] ? "passed" : "failed");
      std::cout << std::endl;
    }
    std::cout << "----------------\n";
  }
  std::cout << name_;
  std::cout << " Unit Tests: " << (test_success ? "passed" : "failed") << std::endl;
}

void TestLookupTable::IntegrationTests() {
  std::cout << name_;
  std::cout << " Integration Tests\n";
}

bool TestLookupTable::UnitTestBobFacsimile(int test_num) {
  bool success = true;
  std::string subtest = "BobFacsimile";

  int ntests = node_[name_][subtest]["test"].size();
  unit_tests_results_[test_num].resize(ntests);

  for (int itest = 0; itest < ntests; ++itest) {
    int ndim = node_[name_][subtest]["test"][itest]["ndim"].as<int>();
    std::vector<double> x[2];
    x[0] = node_[name_][subtest]["test"][itest]["x0"].as<std::vector<double>>();
    x[1] = node_[name_][subtest]["test"][itest]["x1"].as<std::vector<double>>();

    xlh::xlink_params params;
    params.alpha = node_[name_][subtest]["test"][itest]["alpha"].as<double>();
    params.r0 = node_[name_][subtest]["test"][itest]["r0"].as<double>();

    // Run the lookup build
    Init(ndim, x, &xlh::prob_1_2, &params);

    std::vector<std::vector<double>> results;
    for (int i = 0; i < node_[name_][subtest]["test"][itest]["result"].size(); ++i) {
      results.push_back(node_[name_][subtest]["test"][itest]["result"][i].as<std::vector<double>>());
    }
    double tolerance = node_[name_][subtest]["test"][itest]["tolerance"].as<double>();

    unit_tests_results_[test_num][itest] = true;
    for (int i = 0; i < results.size(); ++i) {
      // Check against the table
      for (int j = 0; j < results[i].size(); ++j) {
        double myres = results[i][j];
        double mytest = table_[i * x_[1].size() + j];

        // Test equality
        if (!uth::almost_equal<double>(mytest, myres, tolerance)) {
          std::cout << "Check failed " << i << ", " << j << " " << mytest << " != " << myres << std::endl;
          unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
        }
      }
    }

  }

  return success;
}
