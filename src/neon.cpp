#include "neon.h"

void Neon::Init() {
  MDBead::Init();
  orientation_[1] -=orientation_[1];
}

void NeonSpecies::InitPotentials (system_parameters *params) {
  AddPotential(SID::neon, SID::neon, new LJ126(params->lj_epsilon,params->neon_diameter,space_, 2.5*params->neon_diameter));
  double ar_ne_diameter = 0.5*(params->argon_diameter+params->neon_diameter);
  AddPotential(SID::argon, SID::neon, new LJ126(params->lj_epsilon,ar_ne_diameter,space_, 2.5*ar_ne_diameter));
}

