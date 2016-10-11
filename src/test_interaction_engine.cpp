#include "test_interaction_engine.h"

#include "br_rod.h"
#include "test_helpers.h"
#include "spindle_pole_body.h"

#include <iomanip>
#include <yaml-cpp/yaml.h>

void TestInteractionEngine::InitTestModule(const std::string &filename) {
  name_ = "TestInteractionEngine";
  filename_ = filename;
  std::cout << name_ << " Init Tests\n";

  // Add the tests outlined in the yaml file?
  node_ = YAML::LoadFile(filename_);

  if (node_[name_]["TestInit"]) {
    std::function<bool(int)> f = std::bind(&TestInteractionEngine::TestInit, this, std::placeholders::_1);
    tests_.push_back(f);
    test_names_.push_back("TestInit");
    test_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["TestTetherParticlesMP"]) {
    std::function<bool(int)> f = std::bind(&TestInteractionEngine::TestTetherParticlesMP, this, std::placeholders::_1);
    tests_.push_back(f);
    test_names_.push_back("TestTetherParticlesMP");
    test_results_.push_back(std::vector<bool>());
  }
}

void TestInteractionEngine::RunTests() {
  std::cout << name_ << " Tests ->\n";
  bool test_success = true;

  for (int i = 0; i < tests_.size(); ++i) {
    bool this_success = tests_[i](i);
    test_success &= this_success;

    std::cout << "----------------\n";
    std::cout << "--Test " << test_names_[i] << ": ";
    std::cout << (this_success ? "PASSED" : "FAILED") << std::endl;
    for (int didpass = 0; didpass < test_results_[i].size(); ++didpass) {
      std::cout << "    Test[" << didpass << "] " << (test_results_[i][didpass] ? "passed" : "failed");
      std::cout << std::endl;
    }
    std::cout << "----------------\n";
  }
  std::cout << name_;
  std::cout << " Tests: " << (test_success ? "passed" : "failed") << std::endl;
}

bool TestInteractionEngine::TestInit(int test_num) {
  bool success = true;
  std::string subtest = "TestInit";

  int ntests = (int)node_[name_][subtest]["test"].size(); 

  test_results_[test_num].resize(ntests);

  for (int itest = 0; itest < ntests; ++itest) {
    params_sub_.n_dim = node_[name_][subtest]["test"][itest]["ndim"].as<int>();
    params_sub_.n_periodic = node_[name_][subtest]["test"][itest]["nperiodic"].as<int>();
    space_sub_.Init(&params_sub_, 10);

    Init(space_sub_.GetStruct(), nullptr, nullptr, nullptr, 1.0);

    test_results_[test_num][itest] = true;
    success = true;
  }

  return success;
}

bool TestInteractionEngine::TestTetherParticlesMP(int test_num) {
  bool success = true;
  std::string subtest = "TestTetherParticlesMP";

  int ntests = (int)node_[name_][subtest]["test"].size();
  test_results_[test_num].resize(ntests);

  // Fake out the simples and position map (will need later)
  oid_position_map_ = new std::unordered_map<int, int>();
  simples_ = new std::vector<Simple*>();
  anchors_ = new std::unordered_map<int, std::vector<anchor_t>>();

  for (int itest = 0; itest < ntests; ++itest) {
    simples_->clear();
    oid_position_map_->clear();
    anchors_->clear();

    // Set up the space structure
    params_sub_.n_dim       = node_[name_][subtest]["test"][itest]["ndim"].as<int>();
    params_sub_.n_periodic  = node_[name_][subtest]["test"][itest]["nperiodic"].as<int>();
    space_sub_.Init(&params_sub_, 10);
    space_ = space_sub_.GetStruct();

    // Generate the simples
    // Create a rod and species
    BrRodSpecies *testRodSpecies = new BrRodSpecies();
    testRodSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), anchors_ ,10);
    // Add the rod
    BrRod *testRod = new BrRod(&params_sub_, space_sub_.GetStruct(), 10, SID::br_rod);
    YAML::Node subnode =  node_[name_][subtest]["test"][itest]["rods"][0];
    BrRodSpecies::CreateTestRod(&testRod, ndim_, simples_, oid_position_map_, &subnode);
    testRodSpecies->AddMember(testRod);
    
    // Create an SPB species
    SpindlePoleBodySpecies *testSPBSpecies = new SpindlePoleBodySpecies();
    testSPBSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), anchors_, 10);
    // Add an SPB
    SpindlePoleBody *testSPB = new SpindlePoleBody(&params_sub_, space_sub_.GetStruct(), 10, SID::spb);
    YAML::Node subnode2 = node_[name_][subtest]["test"][itest]["spbs"][0];
    SpindlePoleBodySpecies::CreateTestSPB(&testSPB, ndim_, simples_, oid_position_map_, &subnode2);

    // Check everything quick
    testRod->Dump();
    testSPB->Dump();
    // Print to make sure working
    std::cout << "simples: \n";
    for (int i = 0; i < simples_->size(); ++i) {
      std::cout << "[" << i << "] -> OID " << (*simples_)[i]->GetOID();
      std::cout << " <--> " << (*oid_position_map_)[(*simples_)[i]->GetOID()] << std::endl;
    }

    testRod->SetAnchors(anchors_);
    testSPB->SetAnchors(anchors_);

    // Create a anchor list entry for this spb
    anchors_->insert(std::make_pair(testSPB->GetOID(), std::vector<anchor_t>()));

  }
  return success;
}
