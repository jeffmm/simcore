#ifndef _SIMCORE_ARGON_H_
#define _SIMCORE_ARGON_H_

#include "species.h"
#include "md_bead.h"
#include "auxiliary.h"
#include "lennard_jones_12_6.h"

class Argon : public MDBead {
  public:
    Argon(system_parameters *params, space_struct *space, long seed, SID sid) : MDBead(params, space, seed, sid) {
      diameter_ = params->argon_diameter;
      mass_ = params->argon_mass;
    }
    ~Argon() {}
    Argon(const Argon& that) : MDBead(that) {}
    Argon& operator=(Argon const& that) {MDBead::operator=(that); return *this;} 
    void Init();
    void UpdatePosition() {MDBead::UpdatePosition();}
    void Integrate() {MDBead::Integrate();}
    void UpdateKineticEnergy() {MDBead::UpdateKineticEnergy();}
    double const GetKineticEnergy() {return MDBead::GetKineticEnergy();}
};

class ArgonSpecies : public Species<Argon> {
  protected:
    void InitPotentials (system_parameters *params);

  public:
    ArgonSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::argon);
      InitPotentials(params);
    }
    ~ArgonSpecies() {}
    ArgonSpecies(const ArgonSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }
};

#endif // _SIMCORE_ARGON_H_
