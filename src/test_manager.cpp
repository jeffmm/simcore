#include "test_manager.h"

#include <yaml-cpp/yaml.h>

#define REGISTER_TEST(n,m) test_factory_.register_class<n>(#m);

void TestManager::InitManager(const std::string& filename) {
  YAML::Node node = YAML::LoadFile(filename);

  std::cout << "Welcome to SimCORE Test Manager\n";
  std::cout << "Running tests from ->\n";
  std::cout << " " << filename << std::endl;

  seed_ = node["seed"].as<long>();
}

void TestManager::RegisterTests() {

}

void TestManager::RunTests() {

}
