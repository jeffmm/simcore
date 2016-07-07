// Implementatino for potential manager

#include "potential_manager.h"
#include "helpers.h"
#include "lennard_jones_12_6.h"
#include "wca.h"
#include "sphere_square_well.h"
#include <yaml-cpp/yaml.h>

#define REGISTER_POTENTIAL(n) pot_factory_.register_class<n>(#n);

void
PotentialManager::Init(space_struct *pSpace, char *pFname) {
  space_ = pSpace;
  fname_ = pFname;
  std::cout << "********\n";
  std::cout << "Potentials Load ->\n";
  std::cout << "  file: " << fname_ << std::endl;

  // Register the potentials
  RegisterPotentials();

  ParsePotentials();

  Print();
}

void
PotentialManager::RegisterPotentials() {
  REGISTER_POTENTIAL(LJ126);
  REGISTER_POTENTIAL(WCA);
  REGISTER_POTENTIAL(SphereSquareWell);
}

void
PotentialManager::ParsePotentials() {
  YAML::Node node = YAML::LoadFile(fname_);
  
  npots_ = node["potentials"].size();

  for (int ipot = 0; ipot < npots_; ++ipot) {
    std::string potname = node["potentials"][ipot]["type"].as<std::string>();
    std::string sid1s   = node["potentials"][ipot]["sid1"].as<std::string>();
    std::string sid2s   = node["potentials"][ipot]["sid2"].as<std::string>();
    SID sid1 = StringToSID(sid1s);
    SID sid2 = StringToSID(sid2s);
    // Get the enum type
    PotentialBase* new_pot = (PotentialBase*) pot_factory_.construct(potname);
    new_pot->Init(space_, ipot, node);
    AddPotential(sid1, sid2, new_pot);
  }
}

void PotentialManager::AddPotential(SID sid1, SID sid2, PotentialBase *pot) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (potentials_.count(key1)) return;
  sid_pair key2 = std::make_pair(sid2, sid1);
  if (potentials_.count(key2)) return;
  potentials_[key1] = pot;
}

PotentialBase* PotentialManager::GetPotential(SID sid1, SID sid2) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (potentials_.count(key1)) return potentials_[key1];
  sid_pair key2 = std::make_pair(sid2,sid1);
  if (potentials_.count(key2)) return potentials_[key2];
  return NULL;
}

void PotentialManager::Print() {
  printf("********\n");
  printf("Potentials ->\n");
  for (auto pot=potentials_.begin(); pot!=potentials_.end(); ++pot) {
      printf("{%d,%d} : ",(int) pot->first.first, (int) pot->first.second);
                          pot->second->Print();
  }
}

double PotentialManager::GetMaxRCut() {
  double max_rcut = 0.0;
  for (auto it = potentials_.begin(); it != potentials_.end(); ++it) {
      max_rcut = std::max(max_rcut, it->second->GetRCut());
  }
  return max_rcut;
}

#undef REGISTER_POTENTIAL
