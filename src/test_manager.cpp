#include "test_manager.h"

#include <yaml-cpp/yaml.h>

#include "test_lookup_table.h"
#include "test_xlink_kmc.h"
#include "test_xlink_kmc_moves.h"

#define REGISTER_TEST_MODULE(n) test_module_factory_.register_class<n>(#n);

void TestManager::InitManager(const std::string& filename) {
  YAML::Node node = YAML::LoadFile(filename);
  filename_ = filename;

  RegisterTestModules();

  std::cout << "Welcome to SimCORE Test Manager\n";
  std::cout << "Running tests from ->\n";
  std::cout << " " << filename << std::endl;

  seed_ = node["seed"].as<long>();

  for (auto alltests = test_module_factory_.m_classes.begin(); alltests != test_module_factory_.m_classes.end(); ++alltests) {
    if (node[alltests->first]) {
      std::cout << "TEST MODULE: " << alltests->first << std::endl;
      TestModuleBase *newtest = (TestModuleBase*)test_module_factory_.construct(alltests->first);
      newtest->InitTestModule(node[alltests->first].as<std::string>());
      test_modules_.push_back(newtest);
    }
  }
}

void TestManager::RegisterTestModules() {
  // Register all tests we have access to
  REGISTER_TEST_MODULE(TestLookupTable);
  REGISTER_TEST_MODULE(TestXlinkKMC);
  REGISTER_TEST_MODULE(TestXlinkKMCMoves);
}

void TestManager::RunTestModules() {
  std::cout << "********\n";
  std::cout << "Running test modules\n";
  for (auto testit = test_modules_.begin(); testit != test_modules_.end(); ++testit) {
    //(*testit)->InitTestModule(filename_);
    // Run the unit tests, and report the results
    (*testit)->RunTests();
  }
}
