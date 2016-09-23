// Implementatino for potential manager

#include "potential_manager.h"

#include "boundary_wca.h"
#include "boundary_wca_tip.h"
#include "harmonic.h"
#include "helpers.h"
#include "lennard_jones_12_6.h"
#include "sphere_line_well.h"
#include "sphere_square_well.h"
#include "wca.h"
#include "xlink_harmonic.h"

#include <yaml-cpp/yaml.h>

#define REGISTER_POTENTIAL(n) pot_factory_.register_class<n>(#n);

void
PotentialManager::Init(std::vector<SpeciesBase*> *pSpecies, space_struct *pSpace, char *pFname) {
  space_ = pSpace;
  fname_ = pFname;
  species_ = pSpecies;
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
  REGISTER_POTENTIAL(SphereLineWell);
  REGISTER_POTENTIAL(Harmonic);
  REGISTER_POTENTIAL(XlinkHarmonic);
  REGISTER_POTENTIAL(BoundaryWCA);
  REGISTER_POTENTIAL(BoundaryWCATip);
}

void
PotentialManager::ParsePotentials() {
  YAML::Node node = YAML::LoadFile(fname_);
  
  npots_ = (int)node["potentials"].size();

  for (int ipot = 0; ipot < npots_; ++ipot) {
    std::string potential_type = "external";
    if (node["potentials"][ipot]["type"]) {
      potential_type = node["potentials"][ipot]["type"].as<std::string>();
    }

    if (potential_type.compare("external") == 0) {
      std::string potname = node["potentials"][ipot]["name"].as<std::string>();
      std::cout << potname << std::endl;
      std::string sid1s   = node["potentials"][ipot]["sid1"].as<std::string>();
      std::string sid2s   = node["potentials"][ipot]["sid2"].as<std::string>();
      SID sid1 = StringToSID(sid1s);
      SID sid2 = StringToSID(sid2s);
      // Get the enum type
      PotentialBase* new_pot = (PotentialBase*) pot_factory_.construct(potname);
      new_pot->Init(space_, ipot, node);
      AddPotentialExternal(sid1, sid2, new_pot);
    } else if (potential_type.compare("internal") == 0) {
      std::string potname = node["potentials"][ipot]["name"].as<std::string>();
      std::string sids    = node["potentials"][ipot]["sid"].as<std::string>();
      SID sid = StringToSID(sids); 
      SpeciesBase *sit;
      for (auto msit = species_->begin(); msit != species_->end(); ++msit) {
        if ((*msit)->GetSID() == sid) {
          sit = (*msit);
          break;
        }
      }
      auto internal_pairs = sit->GetInternalPairs();
      PotentialBase *new_pot = (PotentialBase*) pot_factory_.construct(potname);
      new_pot->Init(space_, ipot, node);
      // Add to the internal potential list, which actually owns it
      potential_vec_.push_back(new_pot);
      std::ostringstream potstring;
      potstring << "(" << (int)sid << ", n: " << internal_pairs.size() << ") : ";
      potential_vec_names_types_.push_back(potstring.str());
      for (auto ipair = internal_pairs.begin(); ipair != internal_pairs.end(); ++ipair) {
        AddPotentialInternal(ipair->first, ipair->second, new_pot);
      }
    } else if (potential_type.compare("tether") == 0) {
      // Tethers must be between specific particles
      // Not quite sure how to do this yet, because it might heavily depend on what we're tethering....
      std::cout << "Tethering not quite supported yet, have to figure out use cases, exiting\n";
      exit(1);
    } else if (potential_type.compare("boundary") == 0) {
      // Boundary potentials!
      std::string potname = node["potentials"][ipot]["name"].as<std::string>();
      std::string sids    = node["potentials"][ipot]["sid"].as<std::string>();
      SID sid = StringToSID(sids); 
      PotentialBase *new_pot = (PotentialBase*) pot_factory_.construct(potname);
      new_pot->Init(space_, ipot, node);

      boundaries_[sid] = new_pot;
    } else {
      std::cout << "Potential type " << potential_type << " not yet supported, exiting\n";
      exit(1);
    }
  }
}

void PotentialManager::AddPotentialExternal(SID sid1, SID sid2, PotentialBase *pot) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (potentials_.count(key1)) return;
  sid_pair key2 = std::make_pair(sid2, sid1);
  if (potentials_.count(key2)) return;
  potentials_[key1] = pot;
}

void PotentialManager::AddPotentialInternal(unsigned int oid1, unsigned int oid2, PotentialBase *pot) {
  auto key1 = std::make_pair(oid1, oid2);
  if (internal_potentials_.count(key1)) return;
  auto key2 = std::make_pair(oid2, oid1);
  if (internal_potentials_.count(key2)) return;
  internal_potentials_[key1] = pot;
}

PotentialBase* PotentialManager::GetPotentialExternal(SID sid1, SID sid2) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (potentials_.count(key1)) return potentials_[key1];
  sid_pair key2 = std::make_pair(sid2,sid1);
  if (potentials_.count(key2)) return potentials_[key2];
  return NULL;
}

PotentialBase *PotentialManager::GetPotentialInternal(unsigned int oid1, unsigned int oid2) {
  auto key1 = std::make_pair(oid1, oid2);
  if (internal_potentials_.count(key1)) return internal_potentials_[key1];
  auto key2 = std::make_pair(oid2, oid1);
  if (internal_potentials_.count(key2)) return internal_potentials_[key2];
  return NULL;
}

PotentialBase* PotentialManager::GetPotentialTether(unsigned int oid1, unsigned int oid2) {
  // Get the potential for tethering by OID
  auto key1 = std::make_pair(oid1, oid2);
  if (tethers_.count(key1)) return tethers_[key1];
  auto key2 = std::make_pair(oid2, oid1);
  if (tethers_.count(key2)) return tethers_[key2];
  return NULL;
}

PotentialBase* PotentialManager::GetPotentialBoundary(SID sid) {
  if (boundaries_.count(sid)) return boundaries_[sid];
  return NULL;
}

std::vector<PotentialBase*> PotentialManager::GetAllPotentials() {
  // Needed to get all the potentials for making sure that the kmc
  // potentials (k and requil) are properly synced between configurations
  std::vector<PotentialBase*> allpots;
  for (auto pot = potentials_.begin(); pot != potentials_.end(); ++pot) {
    allpots.push_back(pot->second);
  }
  for (auto pot = internal_potentials_.begin(); pot != internal_potentials_.end(); ++pot) {
    allpots.push_back(pot->second);
  }
  for (auto pot = tethers_.begin(); pot != tethers_.end(); ++pot) {
    allpots.push_back(pot->second);
  }
  return allpots;
}

void PotentialManager::Print() {
  std::cout << "****************\n";
  std::cout << "Potentials ->\n";
  // External potentials
  std::cout << "----------------\n";
  std::cout << "External potentials: \n";
  for (auto pot=potentials_.begin(); pot!=potentials_.end(); ++pot) {
      printf("{%d,%d} : ",(int) pot->first.first, (int) pot->first.second);
                          pot->second->Print();
  }

  // Internal potentials
  std::cout << "----------------\n";
  std::cout << "Internal potentials: \n";
  //for (auto pot=internal_potentials_.begin(); pot != internal_potentials_.end(); ++pot) {
  //  std::cout << "(" << pot->first.first << ", " << pot->first.second << ") : ";
  //  pot->second->Print();
  //}
  for (int ipot = 0; ipot < (int)potential_vec_.size(); ++ipot) {
    std::cout << potential_vec_names_types_[ipot];
    potential_vec_[ipot]->Print();
  }

  // Tethering potentials
  std::cout << "----------------\n";
  std::cout << "Tethering potentials: \n";
  for (auto pot=tethers_.begin(); pot!=tethers_.end(); ++pot) {
    std::cout << "(" << pot->first.first << ", " << pot->first.second << ") : ";
    pot->second->Print();
  }

  // Boundaries
  std::cout << "----------------\n";
  std::cout << "Boundary potentials: \n";
  for (auto pot = boundaries_.begin(); pot != boundaries_.end(); ++pot) {
    std::cout << "(" << (int)pot->first << ", boundary) : ";
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
