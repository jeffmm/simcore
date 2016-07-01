#include "argon.h"

void Argon::Init() {
  MDBead::Init();
}

//void ArgonSpecies::InitPotentials(system_parameters *params) {
  //AddPotential(SID::argon, SID::argon, new LJ126(params->lj_epsilon,params->argon_diameter,space_, 2.5*params->argon_diameter));
  ////AddPotential(SID::argon, SID::argon, new LJ126(params->lj_epsilon,params->argon_diameter,space_, params->argon_rcutoff));
  //double ar_ne_diameter = 0.5*(params->argon_diameter+params->neon_diameter);
  //AddPotential(SID::argon, SID::neon, new LJ126(params->lj_epsilon,ar_ne_diameter,space_, 2.5*ar_ne_diameter));
//}

