#include "test_manager.h"

#include <yaml-cpp/yaml.h>

#include "test_xlink_kmc.h"

#define REGISTER_TEST(n) test_factory_.register_class<n>(#n);

void TestManager::InitManager(const std::string& filename) {
  YAML::Node node = YAML::LoadFile(filename);

  RegisterTests();

  std::cout << "Welcome to SimCORE Test Manager\n";
  std::cout << "Running tests from ->\n";
  std::cout << " " << filename << std::endl;

  seed_ = node["seed"].as<long>();
}

void TestManager::RegisterTests() {
  // Register all tests we have access to
  std::cout << "Registering tests\n";
  REGISTER_TEST(TestXlinkKMC);
}

void TestManager::RunTests() {
  std::cout << "Running tests\n";
  TestBase *xlinktest = (TestBase*)test_factory_.construct("TestXlinkKMC");
  tests_.push_back(xlinktest);

  for (auto testit = tests_.begin(); testit != tests_.end(); ++testit) {
    (*testit)->InitTests();
  }
}
