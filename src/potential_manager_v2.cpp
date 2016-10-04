// Implementatino for potential manager

#include "potential_manager_v2.h"

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

void PotentialManagerV2::Init(std::vector<SpeciesBase*> *pSpecies, space_struct *pSpace, al_set *pAnchors) {
  std::cout << "Potential Manager V2 Init\n";
  space_ = pSpace;
  species_ = pSpecies;
  anchors_ = pAnchors;

  // Register the potentials
  RegisterPotentials();
  std::cout << "Potential Manager V2 Init done\n";
}

void PotentialManagerV2::RegisterPotentials() {
  REGISTER_POTENTIAL(LJ126);
  REGISTER_POTENTIAL(WCA);
  REGISTER_POTENTIAL(SphereSquareWell);
  REGISTER_POTENTIAL(SphereLineWell);
  REGISTER_POTENTIAL(Harmonic);
  REGISTER_POTENTIAL(XlinkHarmonic);
  REGISTER_POTENTIAL(BoundaryWCA);
  REGISTER_POTENTIAL(BoundaryWCATip);
}

// Create the potential, and return the index into the potentials array
int PotentialManagerV2::AddPotential(YAML::Node *subnode) {
  YAML::Node node = *subnode;
  std::string potname = node["name"].as<std::string>();
  PotentialBase *new_pot = (PotentialBase*) pot_factory_.construct(potname);
  new_pot->Init(space_, &node);
  potentials_.push_back(new_pot);
  return (int)potentials_.size() -1;
}

// Print out the information
void PotentialManagerV2::Print() {
  std::cout << "********\n";
  std::cout << "Potentials ->\n";
  // Print out the idx and then each potentials
  for (int ipot = 0; ipot < (int)potentials_.size(); ++ipot) {
    std::cout << "[" << ipot << "] : ";
    potentials_[ipot]->Print();
  }
}
