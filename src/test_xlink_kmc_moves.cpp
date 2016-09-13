#include "test_xlink_kmc_moves.h"

#include "br_rod.h"
#include "test_helpers.h"
#include "xlink.h"
#include "xlink_helpers.h"

#include <iomanip>
#include <yaml-cpp/yaml.h>

void TestXlinkKMCMoves::InitTestModule(const std::string& filename) {
  name_ = "TestXlinkKMCMoves";
  filename_ = filename;
  std::cout << name_ << " Init Tests\n";

  node_ = YAML::LoadFile(filename_);

  if (node_[name_]["KMC_0_1"]) {
    std::function<bool(int)> f = std::bind(&TestXlinkKMCMoves::TestKMC_0_1, this, std::placeholders::_1);
    tests_.push_back(f);
    test_names_.push_back("KMC_0_1");
    test_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["KMC_1_0"]) {
    std::function<bool(int)> f = std::bind(&TestXlinkKMCMoves::TestKMC_1_0, this, std::placeholders::_1);
    tests_.push_back(f);
    test_names_.push_back("KMC_1_0");
    test_results_.push_back(std::vector<bool>());
  }

  sid1_ = SID::xlink;
  sid2_ = SID::br_rod;
}

void TestXlinkKMCMoves::RunTests() {
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

// Tests themselves
bool TestXlinkKMCMoves::TestKMC_0_1(int test_num) {
  bool success = true;
  std::string subtest = "KMC_0_1";

  // Lots of things to do...
  int ntests = (int)node_[name_][subtest]["test"].size();
  test_results_[test_num].resize(ntests);

  // Fake out the simples and position map (will need later)
  oid_position_map_ = new std::unordered_map<int, int>();
  simples_ = new std::vector<Simple*>();
  // Only ever have 1 xlink, 2 heads
  neighbors_ = new nl_list[2];

  for (int itest = 0; itest < ntests; ++itest) {
    ndim_ = node_[name_][subtest]["test"][itest]["ndim"].as<int>();
    params_sub_.n_dim = ndim_;
    params_sub_.delta = node_[name_][subtest]["test"][itest]["delta"].as<double>();
    space_sub_.Init(&params_sub_, 10);
    space_ = space_sub_.GetStruct();

    simples_->clear();
    oid_position_map_->clear();
    neighbors_[0].clear();
    neighbors_[1].clear();

    // Create an xlink and species...
    // Species first
    XlinkSpecies *testXlinkSpecies = new XlinkSpecies();
    testXlinkSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), 10);
    spec1_ = testXlinkSpecies;
    Xlink *testXlink = new Xlink(&params_sub_, space_sub_.GetStruct(), 10, SID::xlink);
    XlinkSpecies::CreateTestXlink(&testXlink, ndim_, simples_, oid_position_map_, filename_, name_, subtest, "xlink", itest);

    testXlinkSpecies->AddMember(testXlink);

    // Set up anything else
    eps_eff_0_1_[0] = 0.5 * node_[name_][subtest]["test"][itest]["concentration_0_1"].as<double>();
    eps_eff_0_1_[1] = eps_eff_0_1_[0];
    alpha_ = node_[name_][subtest]["test"][itest]["alpha"].as<double>();
    on_rate_0_1_[0] = node_[name_][subtest]["test"][itest]["on_rate_0_1"].as<double>();
    on_rate_0_1_[1] = on_rate_0_1_[0];
    rcutoff_0_1_ = 1.0;

    // Neighbors, ugh, also however many rods we happen to have...
    // Loop over the number of them, and create rods with the appropriate attributes
    int nrods = (int)node_[name_][subtest]["test"][itest]["neighbors"].size();
    BrRodSpecies *testRodSpecies = new BrRodSpecies();
    testRodSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), 10);
    for (int irod = 0; irod < nrods; ++irod) {
      // Add the rod
      BrRod *testRod = new BrRod(&params_sub_, space_sub_.GetStruct(), 10, SID::br_rod);
      YAML::Node subnode =  node_[name_][subtest]["test"][itest]["neighbors"][irod];
      BrRodSpecies::CreateTestRod(&testRod, ndim_, simples_, oid_position_map_, &subnode);
      testRodSpecies->AddMember(testRod);

      // Add the neighbors
      neighbor_t new_neighbor;
      new_neighbor.idx_ = (*oid_position_map_)[testRod->GetSimples()[0]->GetOID()];
      new_neighbor.kmc_ = node_[name_][subtest]["test"][itest]["neighbors"][irod]["kmc"].as<double>();
      neighbors_[0].push_back(new_neighbor);
      neighbors_[1].push_back(new_neighbor);
    }

    //Check to make sure everything is working....
    std::cout << "DEBUG check: \n";
    std::cout << "Simples <--> OID:\n";
    for (int i = 0; i < simples_->size(); ++i) {
      std::cout << "[" << i << "] -> OID " << (*simples_)[i]->GetOID();
      std::cout << " <--> " << (*oid_position_map_)[(*simples_)[i]->GetOID()] << std::endl;
    }
    std::cout << "Neighbor list:\n";
    int ineighb = 0;
    for (auto nldx = neighbors_[0].begin(); nldx != neighbors_[0].end(); ++nldx, ++ineighb) {
      std::cout << "[0] = [1] -> [" << ineighb << "]: " << nldx->idx_ << ", kmc: " << std::setprecision(16) << nldx->kmc_ << std::endl;
    }

    // Dump it all
    testXlinkSpecies->Dump();
    testXlinkSpecies->DumpKMC();
    testRodSpecies->Dump();

    KMC_0_1();

    // Test the results
    test_results_[test_num][itest] = true;
    XlinkHead *freehead, *boundhead;
    auto isbound = testXlink->GetBoundHeads(&freehead, &boundhead);

    int result_bound = node_[name_][subtest]["test"][itest]["resulthead"].as<int>();
    double tolerance = node_[name_][subtest]["test"][itest]["tolerance"].as<double>();
    if (isbound.first && !(result_bound == 0)) {
      std::cout << "Bound the wrong head....\n";
      test_results_[test_num][itest] = false && test_results_[test_num][itest];
    }

    double result = node_[name_][subtest]["test"][itest]["result"].as<double>();
    if (!uth::almost_equal<double>(result, boundhead->GetAttach().second, tolerance)) {
      std::cout << "Wrong position, expected: " << std::setprecision(16) << result << ", got: " << boundhead->GetAttach().second << std::endl;
      test_results_[test_num][itest] = false && test_results_[test_num][itest];
    }

    // Clean up after ourselves, yuck
    delete testXlinkSpecies;
    delete testRodSpecies;
  }

  // Clean up
  simples_->clear();
  delete simples_;
  oid_position_map_->clear();
  delete oid_position_map_;

  for (auto suc = test_results_[test_num].begin(); suc != test_results_[test_num].end(); ++suc) {
    success = success && (*suc);
  }
  return success;
}

bool TestXlinkKMCMoves::TestKMC_1_0(int test_num) {
  bool success = true;
  std::string subtest = "KMC_1_0";

  // Lots of things to do...
  int ntests = (int)node_[name_][subtest]["test"].size();
  test_results_[test_num].resize(ntests);

  // Fake out the simples and position map (will need later)
  oid_position_map_ = new std::unordered_map<int, int>();
  simples_ = new std::vector<Simple*>();

  // Run the tests
  for (int itest = 0; itest < ntests; ++itest) {
    ndim_ = node_[name_][subtest]["test"][itest]["ndim"].as<int>();
    params_sub_.n_dim = ndim_;
    params_sub_.delta = node_[name_][subtest]["test"][itest]["delta"].as<double>();
    space_sub_.Init(&params_sub_, 10);
    space_ = space_sub_.GetStruct();

    simples_->clear();
    oid_position_map_->clear();

    // Create a rod and species
    BrRodSpecies *testRodSpecies = new BrRodSpecies();
    testRodSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), 10);
    // Add the rod
    BrRod *testRod = new BrRod(&params_sub_, space_sub_.GetStruct(), 10, SID::br_rod);
    YAML::Node subnode =  node_[name_][subtest]["test"][itest]["rods"][0];
    BrRodSpecies::CreateTestRod(&testRod, ndim_, simples_, oid_position_map_, &subnode);
    testRodSpecies->AddMember(testRod);

    // Create the xlink
    XlinkSpecies *testXlinkSpecies = new XlinkSpecies();
    testXlinkSpecies->InitConfig(&params_sub_, space_sub_.GetStruct(), 10);
    spec1_ = testXlinkSpecies;
    int ncopies = node_[name_][subtest]["test"][itest]["copyfactor_xlink"].as<int>();
    int ninterest = node_[name_][subtest]["test"][itest]["interest_xlink"].as<int>();
    XlinkHead* finalhead;
    for (int ilink = 0; ilink < ncopies; ++ilink) {
      Xlink *testXlink = new Xlink(&params_sub_, space_sub_.GetStruct(), 10, SID::xlink);
      XlinkSpecies::CreateTestXlink(&testXlink, ndim_, simples_, oid_position_map_, filename_, name_, subtest, "xlink", itest, testRod->GetSimples()[0]->GetOID());
      testXlink->UpdateStagePosition(testRod->GetSimples()[0]->GetRigidPosition(),
                                     testRod->GetSimples()[0]->GetRigidOrientation(),
                                     testRod->GetSimples()[0]->GetRigidLength(),
                                     testRod->GetSimples()[0]->GetOID(),
                                     nullptr,
                                     nullptr,
                                     0.0,
                                     0);
      testXlinkSpecies->AddMember(testXlink);
      XlinkHead *freehead, *boundhead;
      auto testbound = testXlink->GetBoundHeads(&freehead, &boundhead);
      std::string rngfile_head = node_[name_][subtest]["test"][itest]["rngfile_head"].as<std::string>();
      freehead->SetRNGState(rngfile_head);
      boundhead->SetRNGState(rngfile_head);
      if (ilink == ninterest) {
        finalhead = boundhead;
      }
    }

    // Load our rng file
    rng_.init(10);
    std::string rngfile = node_[name_][subtest]["test"][itest]["rngfile_kmc"].as<std::string>();
    SetRNGState(rngfile);

    // Load everything else
    on_rate_0_1_[0] = node_[name_][subtest]["test"][itest]["on_rate_0_1"].as<double>();
    on_rate_0_1_[1] = on_rate_0_1_[0];
    rcutoff_0_1_ = 1.0;
    int nbound1[2] = {0, 0};
    nbound1[0] = node_[name_][subtest]["test"][itest]["nbound1"][0].as<int>();
    nbound1[1] = node_[name_][subtest]["test"][itest]["nbound1"][1].as<int>();
    testXlinkSpecies->SetNBound1(nbound1[0], nbound1[1]);
    alpha_ = node_[name_][subtest]["test"][itest]["alpha"].as<double>();

    KMC_1_0();

    finalhead->Dump();
    test_results_[test_num][itest] = true;

    // Get results
    double tolerance = node_[name_][subtest]["test"][itest]["tolerance"].as<double>();
    for (int idim = 0; idim < ndim_; ++idim) {
      double result = node_[name_][subtest]["test"][itest]["result"][idim].as<double>();
      if (!uth::almost_equal<double>(result, finalhead->GetPosition()[idim], tolerance)) {
        std::cout << std::setprecision(16) << "expected: " << result << ", got: " << finalhead->GetPosition()[idim] << std::endl;
        test_results_[test_num][itest] = false && test_results_[test_num][itest];
      }
    }
  }

  // Clean up
  simples_->clear();
  delete simples_;
  oid_position_map_->clear();
  delete oid_position_map_;

  for (auto suc = test_results_[test_num].begin(); suc != test_results_[test_num].end(); ++suc) {
    success = success && (*suc);
  }
  return success;
}
