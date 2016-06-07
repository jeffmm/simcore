#include "neon.h"

// Initialize position, etc
void Neon::Init() {
  MDBead::Init();
  // Orientation is fixed: sets particle color
  orientation_[1] -=orientation_[1];
}

// Initialize interaction potentials for particle
void NeonSpecies::InitPotentials (system_parameters *params) {
  double ne_d = params->neon_diameter;
  double ar_ne_d = 0.5*(ne_d+params->argon_diameter);
  // Neon-neon LJ potential
  AddPotential(SID::neon, SID::neon, 
               new LJ126(params->lj_epsilon,ne_d,space_,2.5*ne_d));
  // Neon-Argon LJ potential
  AddPotential(SID::argon, SID::neon, 
      new LJ126(params->lj_epsilon,ar_ne_d,space_, 2.5*ar_ne_d));
}

