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
  if (node_[name_]["CalcCutoff"]) {
    std::function<bool(void)> f_calc_cutoff = std::bind(&TestXlinkKMC::UnitTestCalcCutoff, this);
    unit_tests_.push_back(f_calc_cutoff);
  }
  if (node_[name_]["Update_0_1"]) {
    std::function<bool(void)> f_update_0_1 = std::bind(&TestXlinkKMC::UnitTestUpdate_0_1, this);
    unit_tests_.push_back(f_update_0_1);
  }
  if (node_[name_]["Update_1_2"]) {
    std::function<bool(void)> f_update_1_2 = std::bind(&TestXlinkKMC::UnitTestUpdate_1_2, this);
    unit_tests_.push_back(f_update_1_2);
  }

  // Create a params, space, etc, and xlink species to test with
  space_sub_.Init(&params_sub_, 10);
  sid2_ = SID::br_rod;
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
    test_success &= (*fit)();
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

  std::vector<bool> sub_success;
  sub_success.reserve(results.size());

  for (int i = 0; i < results.size(); ++i) {
    barrier_weight_ = barrier_weights[i];
    r_equil_ = equilibrium_lengths[i];

    CalcCutoff();

    double res = rcutoff_1_2_;

    if (uth::almost_equal<double>(res, results[i], 1e-8)) {
      std::cout << "    Test[" << i << "] passed\n";
      sub_success[i] = true;
    } else {
      std::cout << "    Test[" << i << "] failed\n";
      sub_success[i] = false;
      success = false;
    }
  }

  return success;
}

bool TestXlinkKMC::UnitTestUpdate_0_1() {
  std::cout << "  Unit Testing XlinkKMC Update_0_1\n";
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

  std::vector<int> n_neighbs;
  std::vector<std::vector<double>> neighbs_vals;
  std::vector<double> results;

  for (int i = 0; i < node_[name_]["Update_0_1"]["results"].size(); ++i) {
    int icount = 0;
    n_neighbs.push_back(node_[name_]["Update_0_1"]["results"][i][icount++].as<double>());
    std::vector<double> values;
    for (int j = 0; j < n_neighbs[i]; ++j) {
      values.push_back(node_[name_]["Update_0_1"]["results"][i][icount++].as<double>());
    }
    neighbs_vals.push_back(values);
    results.push_back(node_[name_]["Update_0_1"]["results"][i][icount++].as<double>());
  }

  std::vector<bool> sub_success;
  sub_success.reserve(results.size());

  // Run the different tests
  for (int i = 0; i < results.size(); ++i) {
    neighbors_[0].clear();
    neighbors_[0].resize(n_neighbs[i]);
    int icount1 = 0;
    for (auto nldx = neighbors_[0].begin(); nldx != neighbors_[0].end(); ++nldx) {
      nldx->kmc_ = neighbs_vals[i][icount1++];
    }
    neighbors_[1].clear();
    neighbors_[1].resize(n_neighbs[i]);
    int icount2 = 0;
    for (auto nldx = neighbors_[1].begin(); nldx != neighbors_[1].end(); ++nldx) {
      nldx->kmc_ = neighbs_vals[i][icount2++];
    }

    // Call Update_0_1
    Update_0_1(testXlink);

    double res = testXlink->GetNExp_0_1();
    if (uth::almost_equal<double>(res, results[i], 1e-8)) {
      std::cout << "    Test[" << i << "] passed\n";
      sub_success[i] = true;
    } else {
      std::cout << "    Test[" << i << "] failed\n";
      sub_success[i] = false;
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

bool TestXlinkKMC::UnitTestUpdate_1_2() {
  std::cout << "  Unit Testing XlinkKMC Update_1_2\n";
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
  neighbors_ = new nl_list[3];
  simples_ = new std::vector<Simple*>();

  eps_eff_1_2_[0] = eps_eff_1_2_[1] = node_[name_]["Update_1_2"]["concentration_1_2"].as<double>();
  on_rate_1_2_[0] = on_rate_1_2_[1] = node_[name_]["Update_1_2"]["on_rate_1_2"].as<double>();
  max_length_ = node_[name_]["Update_1_2"]["max_length"].as<double>();

  // Use the configurators to initialize the xlink and rod
  int ntests = (int)node_[name_]["Update_1_2"]["test"].size();
  std::cout << "   Found " << ntests << " tests\n";

  std::vector<bool> sub_success;
  sub_success.reserve(ntests);

  for (int i = 0; i < ntests; ++i) {
    // Clear the oid map, and the simples
    // Add the base rod
    oid_position_map_->clear();
    simples_->clear();
    double rx_base[3], ux_base[3];
    for (int idim = 0; idim < ndim_; ++idim) {
      rx_base[idim] = node_[name_]["Update_1_2"]["test"][i]["x_rod_base"][idim].as<double>();
      ux_base[idim] = node_[name_]["Update_1_2"]["test"][i]["u_rod_base"][idim].as<double>();
    }
    double l_base = node_[name_]["Update_1_2"]["test"][i]["length_base"].as<double>();
    testRodBase ->InitConfigurator(rx_base, ux_base, l_base);
    //testRodBase->Dump();
    simples_->push_back(testRodBase->GetSimples()[0]);
    (*oid_position_map_)[testRodBase->GetSimples()[0]->GetOID()] = (int)simples_->size() - 1;

    // Add the 'free' rod
    double rx_free[3], ux_free[3];
    for (int idim = 0; idim < ndim_; ++idim) {
      rx_free[idim] = node_[name_]["Update_1_2"]["test"][i]["x_rod_free"][idim].as<double>();
      ux_free[idim] = node_[name_]["Update_1_2"]["test"][i]["u_rod_free"][idim].as<double>();
    }
    double l_free = node_[name_]["Update_1_2"]["test"][i]["length_free"].as<double>();
    testRodFree->InitConfigurator(rx_free, ux_free, l_free);
    //testRodFree->Dump();
    simples_->push_back(testRodFree->GetSimples()[0]);
    (*oid_position_map_)[testRodFree->GetSimples()[0]->GetOID()] = (int)simples_->size() - 1;

    // Add the xlink
    double xx[3] = {0.0, 0.0, 0.0}; // we will set based upon the attached rod
    testXlink->InitConfigurator(xx, 1.0);
    int head = node_[name_]["Update_1_2"]["test"][i]["head"].as<double>();
    double crosspos = node_[name_]["Update_1_2"]["test"][i]["cross_position"].as<double>();
    auto bindhead = testXlink->GetHeads()->begin() + head;
    bindhead->Attach(testRodBase->GetOID(), crosspos);
    bindhead->SetBound(true);
    testXlink->CheckBoundState();

    // Set the location
    for (int i = 0; i < ndim_; ++i) {
      xx[i] = rx_base[i] - 0.5 * ux_base[i] * l_base + crosspos * ux_base[i];
    }

    XlinkHead *freehead, *boundhead;
    auto isbound = testXlink->GetBoundHeads(&freehead, &boundhead);
    freehead->SetPosition(xx);
    boundhead->SetPosition(xx);
    testXlink->SetPosition(xx);
    simples_->push_back(freehead);
    (*oid_position_map_)[freehead->GetOID()] = (int)simples_->size() - 1;
    //testXlink->Dump();
    //testXlink->DumpKMC();

    // Just add the neighbor list by hand (ugh)
    auto freehead_idx = (*oid_position_map_)[freehead->GetOID()];
    neighbors_[freehead_idx].clear();
    neighbors_[freehead_idx].resize(1);
    for (auto nldx = neighbors_[freehead_idx].begin(); nldx != neighbors_[freehead_idx].end(); ++nldx) {
      nldx->idx_ = (*oid_position_map_)[testRodFree->GetSimples()[0]->GetOID()];
    }

    //print the neighbor list
    /*std::cout << "neighbor list: \n";
    for (int i = 0; i < 2; ++i) {
      int iset = -1;
      for (auto nldx = neighbors_[i].begin(); nldx != neighbors_[i].end(); ++nldx) {
        iset++;
        printf("[%d] -> [%d] = %d\n", i, iset, nldx->idx_);
      }
    }*/

    // print a bunch of info to make sure it's right
    // oid map
    /*std::cout << "oid position map: \n";
    for (auto it = (*oid_position_map_).begin(); it != (*oid_position_map_).end(); it++) {
      printf("%d <--> %d\n", it->first, it->second);
    }
    std::cout << "simples:\n";
    for (int i = 0; i < simples_->size(); ++i) {
      printf("[%d] -> oid(%d)\n", i, (*simples_)[i]->GetOID());
    }*/

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
      std::cout << "    Test[" << i << "] passed -> ";
      std::cout << "expected: " << this_result << ", got: " << res << " (tol: ";
      std::cout << tolerance << ")\n";
      sub_success[i] = true;
    } else {
      std::cout << "    Test[" << i << "] failed -> ";
      std::cout << "expected: " << this_result << ", got: " << res << std::endl;
      sub_success[i] = false;
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
