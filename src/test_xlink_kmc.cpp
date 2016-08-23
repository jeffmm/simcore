#include "test_xlink_kmc.h"

#include "br_rod.h"
#include "test_helpers.h"
#include "xlink.h"

#include <iomanip>
#include <yaml-cpp/yaml.h>

void TestXlinkKMC::InitTestModule(const std::string& filename) {
  name_ = "TestXlinkKMC";
  std::cout << name_ << " Init Tests\n";

  // Add the tests outlined in the yaml file?
  node_ = YAML::LoadFile(filename);

  // Test out functional
  if (node_[name_]["FailTest"]) {
    std::function<bool(int)> f_failtest = std::bind(&TestXlinkKMC::FailTest, this, std::placeholders::_1);
    unit_tests_.push_back(f_failtest);
    unit_tests_names_.push_back("FailTest");
    unit_tests_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["CalcCutoff"]) {
    std::function<bool(int)> f_calc_cutoff = std::bind(&TestXlinkKMC::UnitTestCalcCutoff, this, std::placeholders::_1);
    unit_tests_.push_back(f_calc_cutoff);
    unit_tests_names_.push_back("CalcCutoff");
    unit_tests_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["Update_0_1"]) {
    std::function<bool(int)> f_update_0_1 = std::bind(&TestXlinkKMC::UnitTestUpdate_0_1, this, std::placeholders::_1);
    unit_tests_.push_back(f_update_0_1);
    unit_tests_names_.push_back("Update_0_1");
    unit_tests_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["Update_1_2"]) {
    std::function<bool(int)> f_update_1_2 = std::bind(&TestXlinkKMC::UnitTestUpdate_1_2, this, std::placeholders::_1);
    unit_tests_.push_back(f_update_1_2);
    unit_tests_names_.push_back("Update_1_2");
    unit_tests_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["Detach_1_0"]) {
    std::function<bool(int)> f_detach_1_0 = std::bind(&TestXlinkKMC::UnitTestDetach_1_0, this, std::placeholders::_1);
    unit_tests_.push_back(f_detach_1_0);
    unit_tests_names_.push_back("Detach_1_0");
    unit_tests_results_.push_back(std::vector<bool>());
  }
  if (node_[name_]["Detach_2_1"]) {
    std::function<bool(int)> f_detach_2_1 = std::bind(&TestXlinkKMC::UnitTestDetach_2_1, this, std::placeholders::_1);
    unit_tests_.push_back(f_detach_2_1);
    unit_tests_names_.push_back("Detach_2_1");
    unit_tests_results_.push_back(std::vector<bool>());
  }

  // Create a params, space, etc, and xlink species to test with
  space_sub_.Init(&params_sub_, 10);
  sid1_ = SID::xlink;
  sid2_ = SID::br_rod;
}

void TestXlinkKMC::RunTests() {
  std::cout << name_ << " Run Tests\n";
  UnitTests();
  IntegrationTests();
}

bool TestXlinkKMC::FailTest(int test_num) {
  // Basic test to make sure that failing works properly
  int ntests = 2;
  bool success = true;

  unit_tests_results_[test_num].resize(ntests);
  
  for (int i = 0; i < ntests; ++i) {
    if (i == 0) {
      unit_tests_results_[test_num][i] = true;
    } else {
      unit_tests_results_[test_num][i] = false;
      success = false;
    }
  }

  return success;
}

void TestXlinkKMC::UnitTests() {
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

void TestXlinkKMC::IntegrationTests() {
  std::cout << name_;
  std::cout << " Integration Tests\n";
}

bool TestXlinkKMC::UnitTestCalcCutoff(int test_num) {
  bool success = true;

  // Grab the YAML file information
  eps_eff_1_2_[0] = eps_eff_1_2_[1] = node_[name_]["CalcCutoff"]["concentration_1_2"].as<double>();
  k_stretch_  = node_[name_]["CalcCutoff"]["spring_constant"].as<double>();
  max_length_ = node_[name_]["CalcCutoff"]["max_length"].as<double>();

  int ntests = node_[name_]["CalcCutoff"]["test"].size();

  unit_tests_results_[test_num].resize(ntests);

  for (int itest = 0; itest < ntests; ++itest) {
    barrier_weight_ = node_[name_]["CalcCutoff"]["test"][itest]["barrier_weight"].as<double>();
    r_equil_ = node_[name_]["CalcCutoff"]["test"][itest]["r_equil"].as<double>();

    CalcCutoff();

    double res = rcutoff_1_2_;

    double result = node_[name_]["CalcCutoff"]["test"][itest]["result"].as<double>();
    double tolerance = node_[name_]["CalcCutoff"]["test"][itest]["tolerance"].as<double>();

    if (uth::almost_equal<double>(res, result, tolerance)) {
      unit_tests_results_[test_num][itest] = true;
    } else {
      unit_tests_results_[test_num][itest] = false;
      success = false;
    }
  }

  return success;
}

bool TestXlinkKMC::UnitTestUpdate_0_1(int test_num) {
  bool success = true;

  // Lots of things to set up for the update test
  params_sub_.xlink_diameter = 1.0;
  params_sub_.delta = 0.00025;

  // Set up the various quantities needed
  Xlink *testXlink = new Xlink(&params_sub_, space_sub_.GetStruct(), 10, SID::xlink);
  oid_position_map_ = new std::unordered_map<int, int>();
  neighbors_ = new nl_list[2];

  // Set the oid map
  auto xheads = testXlink->GetHeads();
  auto head0 = &(*xheads)[0];
  auto head1 = &(*xheads)[1];
  (*oid_position_map_)[head0->GetOID()] = 0;
  (*oid_position_map_)[head1->GetOID()] = 1;

  // Get values
  eps_eff_0_1_[0] = eps_eff_0_1_[1] = 0.5 * node_[name_]["Update_0_1"]["concentration_0_1"].as<double>();
  on_rate_0_1_[0] = on_rate_0_1_[1] = node_[name_]["Update_0_1"]["on_rate_0_1"].as<double>();
  alpha_ = node_[name_]["Update_0_1"]["alpha"].as<double>();

  int ntests = node_[name_]["Update_0_1"]["test"].size();
  unit_tests_results_[test_num].resize(ntests);

  // Run the different tests
  for (int i = 0; i < ntests; ++i) {
    int n_neighbors = node_[name_]["Update_0_1"]["test"][i]["n_neighbors"].as<int>();
    neighbors_[0].clear();
    neighbors_[0].resize(n_neighbors);
    int ineighb = 0;
    for (auto nldx = neighbors_[0].begin(); nldx != neighbors_[0].end(); ++nldx) {
      nldx->kmc_ = node_[name_]["Update_0_1"]["test"][i]["values"][ineighb].as<double>();
      ineighb++;
    }
    neighbors_[1].clear();
    neighbors_[1].resize(n_neighbors);
    ineighb = 0;
    for (auto nldx = neighbors_[1].begin(); nldx != neighbors_[1].end(); ++nldx) {
      nldx->kmc_ = node_[name_]["Update_0_1"]["test"][i]["values"][ineighb].as<double>();
      ineighb++;
    }

    // Call Update_0_1
    Update_0_1(testXlink);

    double res = testXlink->GetNExp_0_1();

    double result = node_[name_]["Update_0_1"]["test"][i]["result"].as<double>();
    double tolerance = node_[name_]["Update_0_1"]["test"][i]["tolerance"].as<double>();
    if (uth::almost_equal<double>(res, result, tolerance)) {
      unit_tests_results_[test_num][i] = true;
    } else {
      unit_tests_results_[test_num][i] = false;
      success = false;
    }
  }

  // Cleanup
  if (testXlink) {
    delete testXlink;
  }
  if (oid_position_map_) {
    delete oid_position_map_;
  }
  if (neighbors_) {
    delete[] neighbors_;
  }

  return success;
}

void TestXlinkKMC::CreateTestRod(BrRod **rod,
                                 const std::string &unitname,
                                 const std::string &rodname,
                                 int itest) {
  BrRod *mrod = *rod;
  // Load the rod
  double xr[3] = {0.0, 0.0, 0.0};
  double ur[3] = {0.0, 0.0, 0.0};
  std::ostringstream xrodname;
  xrodname << "x_" << rodname;
  std::ostringstream urodname;
  urodname << "u_" << rodname;
  std::ostringstream lrodname;
  lrodname << "l_" << rodname;
  //std::cout << "Checking location\n";
  for (int idim = 0; idim < ndim_; ++idim) {
    xr[idim] = node_[name_][unitname.c_str()]["test"][itest][xrodname.str().c_str()][idim].as<double>();
    ur[idim] = node_[name_][unitname.c_str()]["test"][itest][urodname.str().c_str()][idim].as<double>();
  }
  //std::cout << "Checking length\n";
  double lrod = node_[name_][unitname.c_str()]["test"][itest][lrodname.str().c_str()].as<double>();
  //std::cout << "TEST ROD CREATE: \n";
  //std::cout << std::setprecision(16) << "x: (" << xr[0] << ", " << xr[1] << ", " << xr[2] << "), ";
  //std::cout << std::setprecision(16) << "u: (" << ur[0] << ", " << ur[1] << ", " << ur[2] << "), ";
  //std::cout << std::setprecision(16) << "l: " << lrod << std::endl;
  mrod->InitConfigurator(xr, ur, lrod);
  //mrod->Dump();

  // Add the simples to the correct location and the oid position map
  std::vector<Simple*> sim_vec = mrod->GetSimples();
  for (int i = 0; i < sim_vec.size(); ++i) {
    simples_->push_back(sim_vec[i]);
    (*oid_position_map_)[sim_vec[i]->GetOID()] = simples_->size() -1;
  }

  // Print to make sure working
  /*std::cout << "simples: \n";
  for (int i = 0; i < simples_->size(); ++i) {
    std::cout << "[" << i << "] -> OID " << (*simples_)[i]->GetOID();
    std::cout << " <--> " << (*oid_position_map_)[(*simples_)[i]->GetOID()] << std::endl;
  }*/
}

void TestXlinkKMC::CreateTestXlink(Xlink **mxit,
                                   const std::string &unitname,
                                   const std::string &xname,
                                   int itest,
                                   int attachoid) {
  Xlink *xit = *mxit;
  double xx[3] = {0.0, 0.0, 0.0};
  xit->InitConfigurator(xx, 1.0);
  //std::cout << "Xlink base name: " << xname << std::endl;
  std::ostringstream headname;
  headname << "head_" << xname;
  //std::cout << "head name: " << headname.str() << std::endl;
  int head = node_[name_][unitname.c_str()]["test"][itest][headname.str().c_str()].as<int>();
  std::ostringstream crossposname;
  crossposname << "crosspos_" << xname;
  //std::cout << "crosspos name: " << crossposname.str() << std::endl;
  double crosspos = node_[name_][unitname.c_str()]["test"][itest][crossposname.str().c_str()].as<double>();

  xit->BindHeadSingle(head, crosspos, attachoid);
  //xit->Dump();
  //xit->DumpKMC();

  // simples and oid stuff
  std::vector<Simple*> sim_vec = xit->GetSimples();
  for (int i = 0; i < sim_vec.size(); ++i) {
    simples_->push_back(sim_vec[i]);
    (*oid_position_map_)[sim_vec[i]->GetOID()] = simples_->size() -1;
  }

  // Print to make sure working
  /*std::cout << "simples: \n";
  for (int i = 0; i < simples_->size(); ++i) {
    std::cout << "[" << i << "] -> OID " << (*simples_)[i]->GetOID();
    std::cout << " <--> " << (*oid_position_map_)[(*simples_)[i]->GetOID()] << std::endl;
  }*/
}

bool TestXlinkKMC::UnitTestUpdate_1_2(int test_num) {
  bool success = true;

  // Lots of things to set up for the update test
  params_sub_.xlink_diameter = 1.0;
  params_sub_.delta = 0.00025;
  params_sub_.n_dim = 3;

  space_sub_.Init(&params_sub_, 10);
  space_ = space_sub_.GetStruct();
  ndim_ = 3;
  nperiodic_ = 0;

  // Set up the various quantities needed
  BrRod *testRodBase = new BrRod(&params_sub_, space_sub_.GetStruct(), 10, SID::br_rod);
  BrRod *testRodFree = new BrRod(&params_sub_, space_sub_.GetStruct(), 10, SID::br_rod);
  Xlink *testXlink = new Xlink(&params_sub_, space_sub_.GetStruct(), 10, SID::xlink);
  oid_position_map_ = new std::unordered_map<int, int>();
  neighbors_ = new nl_list[4];
  simples_ = new std::vector<Simple*>();

  eps_eff_1_2_[0] = eps_eff_1_2_[1] = node_[name_]["Update_1_2"]["concentration_1_2"].as<double>();
  on_rate_1_2_[0] = on_rate_1_2_[1] = node_[name_]["Update_1_2"]["on_rate_1_2"].as<double>();
  max_length_ = node_[name_]["Update_1_2"]["max_length"].as<double>();

  // Use the configurators to initialize the xlink and rod
  int ntests = (int)node_[name_]["Update_1_2"]["test"].size();

  unit_tests_results_[test_num].resize(ntests);

  for (int i = 0; i < ntests; ++i) {
    // Clear the oid map, and the simples
    // Add the base rod
    oid_position_map_->clear();
    simples_->clear();
    CreateTestRod(&testRodBase, "Update_1_2", "rod_base", i);
    CreateTestRod(&testRodFree, "Update_1_2", "rod_free", i);

    // Add the xlink
    CreateTestXlink(&testXlink, "Update_1_2", "xlink", i, testRodBase->GetSimples()[0]->GetOID());
    double xx[3] = {0.0, 0.0, 0.0}; // we will set based upon the attached rod

    // Set the location
    auto rx_base = testRodBase->GetSimples()[0]->GetRigidPosition();
    auto ux_base = testRodBase->GetSimples()[0]->GetRigidOrientation();
    auto l_base = testRodBase->GetSimples()[0]->GetRigidLength();

    XlinkHead *freehead, *boundhead;
    auto isbound = testXlink->GetBoundHeads(&freehead, &boundhead);
    auto crosspos = boundhead->GetAttach().second;
    for (int i = 0; i < ndim_; ++i) {
      xx[i] = rx_base[i] - 0.5 * ux_base[i] * l_base + crosspos * ux_base[i];
    }
    freehead->SetPosition(xx);
    boundhead->SetPosition(xx);
    testXlink->SetPosition(xx);
    //testXlink->Dump();
    //testXlink->DumpKMC();

    // Just add the neighbor list by hand (ugh)
    auto freehead_idx = (*oid_position_map_)[freehead->GetOID()];
    neighbors_[freehead_idx].clear();
    neighbors_[freehead_idx].resize(1);
    for (auto nldx = neighbors_[freehead_idx].begin(); nldx != neighbors_[freehead_idx].end(); ++nldx) {
      nldx->idx_ = (*oid_position_map_)[testRodFree->GetSimples()[0]->GetOID()];
    }

    r_equil_ = node_[name_]["Update_1_2"]["test"][i]["r_equil"].as<double>();
    k_stretch_ = node_[name_]["Update_1_2"]["test"][i]["k_stretch"].as<double>();
    barrier_weight_ = node_[name_]["Update_1_2"]["test"][i]["barrier_weight"].as<double>();

    // Run the update
    CalcCutoff();
    BuildTables();
    Update_1_2(testXlink);

    auto res = testXlink->GetNExp_1_2();

    double this_result = node_[name_]["Update_1_2"]["test"][i]["result"].as<double>();
    double tolerance = node_[name_]["Update_1_2"]["test"][i]["tolerance"].as<double>();
    if (uth::almost_equal<double>(res, this_result, tolerance)) {
      unit_tests_results_[test_num][i] = true;
    } else {
      unit_tests_results_[test_num][i] = false;
      success = false;
    }

    // Some interstage cleanup
    freehead->SetBound(false);
    freehead->Attach(-1, 0.0);
    boundhead->SetBound(false);
    boundhead->Attach(-1, 0.0);
    testXlink->CheckBoundState();
  }

  // Cleanup
  if (testXlink) {
    delete testXlink;
  }
  if (testRodBase) {
    delete testRodBase;
  }
  if (testRodFree) {
    delete testRodFree;
  }
  if (oid_position_map_) {
    delete oid_position_map_;
  }
  if (neighbors_) {
    delete[] neighbors_;
  }
  if (simples_) {
    simples_->clear();
    delete simples_;
  }

  return success;
}

bool TestXlinkKMC::UnitTestDetach_1_0(int test_num) {
  bool success = true;

  // Use the configurators to initialize the xlink and rod
  int ntests = (int)node_[name_]["Detach_1_0"]["test"].size();

  // Fake out the position map and simples
  oid_position_map_ = new std::unordered_map<int, int>();
  simples_ = new std::vector<Simple*>();

  unit_tests_results_[test_num].resize(ntests);
  
  for (int itest = 0; itest < ntests; ++itest) {
    // Clean up everything from previous tests
    oid_position_map_->clear();
    simples_->clear();
    // Set up the space type structs
    ndim_ = node_[name_]["Detach_1_0"]["test"][itest]["ndim"].as<int>();
    params_sub_.n_dim = ndim_;
    space_sub_.Init(&params_sub_, 10);
    space_ = space_sub_.GetStruct();

    rcutoff_0_1_ = node_[name_]["Detach_1_0"]["test"][itest]["rcutoff_0_1"].as<double>();

    long x_rng = node_[name_]["Detach_1_0"]["test"][itest]["x_rng"].as<long>();

    BrRod *testRod = new BrRod(&params_sub_, space_sub_.GetStruct(), x_rng, SID::br_rod);
    Xlink *testXlink = new Xlink(&params_sub_, space_sub_.GetStruct(), x_rng, SID::xlink);

    CreateTestRod(&testRod, "Detach_1_0", "rod", itest);

    // Load the xlink
    CreateTestXlink(&testXlink, "Detach_1_0", "xlink", itest, testRod->GetSimples()[0]->GetOID());
    double xx[3] = {0.0, 0.0, 0.0};

    XlinkHead *freehead, *boundhead;
    auto isbound = testXlink->GetBoundHeads(&freehead, &boundhead);
    auto xr = testRod->GetSimples()[0]->GetRigidPosition();
    auto ur = testRod->GetSimples()[0]->GetRigidOrientation();
    auto lrod = testRod->GetSimples()[0]->GetRigidLength();
    auto crosspos = boundhead->GetAttach().second;

    // Set the location
    for (int i = 0; i < ndim_; ++i) {
      xx[i] = xr[i] - 0.5 * ur[i] * lrod + crosspos * ur[i];
    }
    freehead->SetPosition(xx);
    boundhead->SetPosition(xx);
    testXlink->SetPosition(xx);
    simples_->push_back(boundhead);
    (*oid_position_map_)[boundhead->GetOID()] = (int)simples_->size() - 1;
    //testXlink->Dump();
    //testXlink->DumpKMC();

    double original_position[3] = {0.0, 0.0, 0.0};
    double orig_results[3] = {0.0, 0.0, 0.0};
    double dr2_orig = 0.0;
    std::copy(boundhead->GetRigidPosition(), boundhead->GetRigidPosition()+ndim_, original_position);
    for (int idim = 0; idim < ndim_; ++idim) {
      orig_results[idim] = node_[name_]["Detach_1_0"]["test"][itest]["results_original_pos"][idim].as<double>();
      dr2_orig += SQR(original_position[idim] - orig_results[idim]);
    }
    // Detach!
    Detach_1_0(testXlink, freehead, boundhead);

    // Check various things
    unit_tests_results_[test_num][itest] = true;
    if (freehead->GetBound()) {
      unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
    }
    if (boundhead->GetBound()) {
      unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
    }
    if (testXlink->GetBoundState() != unbound) {
      unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
    }

    double res[3] = {0.0, 0.0, 0.0};
    double dr_loc = 0.0;
    for (int idim = 0; idim < ndim_; ++idim) {
      res[idim] = node_[name_]["Detach_1_0"]["test"][itest]["results"][idim].as<double>();
      dr_loc += SQR(res[idim] - boundhead->GetRigidPosition()[idim]);
    }

    double tolerance = node_[name_]["Detach_1_0"]["test"][itest]["tolerance"].as<double>();
    if (!uth::almost_equal(dr2_orig, 0.0, tolerance)) {
      unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
    }
    if (!uth::almost_equal(dr_loc, 0.0, tolerance)) {
      unit_tests_results_[test_num][itest] = false && unit_tests_results_[test_num][itest];
    }

    // Cleanup
    if (testRod) {
      delete testRod;
    }
    if (testXlink) {
      delete testXlink;
    }
  }

  if (simples_) {
    delete simples_;
  }
  if (oid_position_map_) {
    delete oid_position_map_;
  }

  for (auto suc = unit_tests_results_[test_num].begin(); suc != unit_tests_results_[test_num].end(); ++suc) {
    success = success && (*suc);
  }
  return success;
}

bool TestXlinkKMC::UnitTestDetach_2_1(int test_num) {
  bool success = true;

  // Use the configurators to initialize the xlink and rod
  int ntests = (int)node_[name_]["Detach_1_0"]["test"].size();

  // Fake out the position map and simples
  oid_position_map_ = new std::unordered_map<int, int>();
  simples_ = new std::vector<Simple*>();

  unit_tests_results_[test_num].resize(ntests);

  for (int itest = 0; itest < ntests; ++itest) {
    // Clean up everything from previous tests
    oid_position_map_->clear();
    simples_->clear();
    // Set up the space type structs
    ndim_ = node_[name_]["Detach_1_0"]["test"][itest]["ndim"].as<int>();
    params_sub_.n_dim = ndim_;
    space_sub_.Init(&params_sub_, 10);
    space_ = space_sub_.GetStruct();
  }

  // Cleanup at end
  if (simples_) {
    delete simples_;
  }
  if (oid_position_map_) {
    delete oid_position_map_;
  }

  return success;
}
