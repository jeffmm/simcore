#include "test_manager.h"

#include <yaml-cpp/yaml.h>

#include "test_xlink_kmc.h"

#define REGISTER_TEST(n) test_factory_.register_class<n>(#n);

void TestManager::InitManager(const std::string& filename) {
  YAML::Node node = YAML::LoadFile(filename);
  filename_ = filename;

  RegisterTests();

  std::cout << "Welcome to SimCORE Test Manager\n";
  std::cout << "Running tests from ->\n";
  std::cout << " " << filename << std::endl;

  seed_ = node["seed"].as<long>();

  for (auto alltests = test_factory_.m_classes.begin(); alltests != test_factory_.m_classes.end(); ++alltests) {
    if (node[alltests->first]) {
      std::cout << "Adding test " << alltests->first << std::endl;
      TestBase *newtest = (TestBase*)test_factory_.construct(alltests->first);
      tests_.push_back(newtest);
    }
  }
}

void TestManager::RegisterTests() {
  // Register all tests we have access to
  REGISTER_TEST(TestXlinkKMC);
}

void TestManager::RunTests() {
  std::cout << "Running tests\n";
  for (auto testit = tests_.begin(); testit != tests_.end(); ++testit) {
    (*testit)->InitTests(filename_);
    // Run the unit tests, and report the results
    (*testit)->UnitTests();
  }
}
