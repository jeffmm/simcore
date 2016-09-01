// Implementatino for potential manager

#include "potential_manager.h"

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
}

void
PotentialManager::ParsePotentials() {
  YAML::Node node = YAML::LoadFile(fname_);
  
  npots_ = (int)node["potentials"].size();

  for (int ipot = 0; ipot < npots_; ++ipot) {
    bool is_tether = false;
    if (node["potentials"][ipot]["tether"]) {
      is_tether = node["potentials"][ipot]["tether"].as<bool>();
    }

    if (is_tether) {
      std::string potname = node["potentials"][ipot]["type"].as<std::string>();
      std::string tether_type = node["potentials"][ipot]["tether_type"].as<std::string>();

      // Only do xlink for now
      if (tether_type.compare("xlink") == 0) {
        // Create the tethering.  Need the object of the members of species to give us
        // back the OIDs of the pairs that we're connecting
        std::string sids = node["potentials"][ipot]["sid"].as<std::string>();
        SID sid = StringToSID(sids); 
        SpeciesBase *sit;
        // Find the species
        for (auto msit = species_->begin(); msit != species_->end(); ++msit) {
          if ((*msit)->GetSID() == sid) {
            sit = (*msit);
            break;
          }
        }
        auto internal_pairs = sit->GetInternalPairs();
        PotentialBase* new_pot = (PotentialBase*) pot_factory_.construct(potname);
        new_pot->Init(space_, ipot, node);
        // Add to the internal potential list
        internal_potentials_.push_back(new_pot);
        int pot_idx = internal_potentials_.size() - 1;
        for (auto ipair = internal_pairs.begin(); ipair != internal_pairs.end(); ++ipair) {
          tethers_[*ipair] = pot_idx;
        }
      } else {
        std::cout << "Tethering type: " << tether_type << " not yet supported, exiting\n";
        exit(1);
      }
    } else {
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

PotentialBase* PotentialManager::GetPotentialTether(unsigned int oid1, unsigned int oid2) {
  std::cout << "GetPotentialTether[" << oid1 << ", " << oid2 << "]\n";
  // Get the potential for tethering by OID
  auto key1 = std::make_pair(oid1, oid2);
  if (tethers_.count(key1)) return internal_potentials_[tethers_[key1]];
  return NULL;
}

void PotentialManager::Print() {
  std::cout << "****************\n";
  std::cout << "Potentials ->\n";
  std::cout << "----------------\n";
  std::cout << "External potentials: \n";
  for (auto pot=potentials_.begin(); pot!=potentials_.end(); ++pot) {
      printf("{%d,%d} : ",(int) pot->first.first, (int) pot->first.second);
                          pot->second->Print();
  }

  // Internal potentials
  std::cout << "----------------\n";
  std::cout << "Internal potentials: \n";
  for (auto ipair = tethers_.begin(); ipair != tethers_.end(); ++ipair) {
    std::cout << "(" << ipair->first.first << ", " << ipair->first.second << ") : ";
    auto mypot = internal_potentials_[ipair->second];
    mypot->Print();
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
