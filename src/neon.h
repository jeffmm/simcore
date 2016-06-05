#ifndef _CYTOSCORE_NEON_H_
#define _CYTOSCORE_NEON_H_

#include "species.h"
#include "md_bead.h"
#include "auxiliary.h"
#include "test_potential.h"
#include "lennard_jones_12_6.h"

class Neon : public MDBead {
  public:
    Neon(system_parameters *params, space_struct *space, long seed, SID sid) : MDBead(params, space, seed, sid) {
      diameter_ = params->neon_diameter;
      mass_ = params->neon_mass;
    }
    ~Neon() {}
    Neon(const Neon& that) : MDBead(that) {}
    Neon& operator=(Neon const& that) {MDBead::operator=(that); return *this;} 
    void Init();
    void UpdatePosition() {MDBead::UpdatePosition();}
    void Integrate() {MDBead::Integrate();}
    void UpdateKineticEnergy() {MDBead::UpdateKineticEnergy();}
    double const GetKineticEnergy() {return MDBead::GetKineticEnergy();}
};

class NeonSpecies : public Species<Neon> {
  protected:
    void InitPotentials (system_parameters *params);
  public:
    NeonSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::neon);
      InitPotentials(params);
    }
    ~NeonSpecies() {}
    NeonSpecies(const NeonSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

};

#endif // _CYTOSCORE_NEON_H_
