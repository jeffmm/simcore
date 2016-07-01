#ifndef _SIMCORE_NEON_H_
#define _SIMCORE_NEON_H_

#include "species.h"
#include "md_bead.h"
#include "auxiliary.h"
#include "lennard_jones_12_6.h"

class Neon : public MDBead {
  public:
    // Constructors, destructors, copy-constructors
    Neon(system_parameters *params, space_struct *space, long seed, 
        SID sid) : MDBead(params, space, seed, sid) {
      // Set parameters unique to Neon
      diameter_ = params->neon_diameter;
      mass_ = params->neon_mass;
    }
    ~Neon() {}
    Neon(const Neon& that) : MDBead(that) {}
    Neon& operator=(Neon const& that) {
      MDBead::operator=(that); return *this;
    }
    // Define simple particle virtual functions for Neon
    void Init();
    void UpdatePosition() {MDBead::UpdatePosition();}
    void Integrate() {MDBead::Integrate();}
    void UpdateKineticEnergy() {MDBead::UpdateKineticEnergy();}
    double const GetKineticEnergy() {return MDBead::GetKineticEnergy();}
};

class NeonSpecies : public Species<Neon> {
  protected:
    //void InitPotentials (system_parameters *params);
  public:
    NeonSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      // Set species ID for Neon
      SetSID(SID::neon);
      //InitPotentials(params);
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

#endif // _SIMCORE_NEON_H_
