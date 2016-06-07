#ifndef _CYTOSCORE_MD_BEAD_H_
#define _CYTOSCORE_MD_BEAD_H_

#include "species.h"
#include "bead.h"
#include "auxiliary.h"
#include "lennard_jones_12_6.h"

class MDBead : public Bead {
  protected:
    double mass_;
  public:
    MDBead(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Bead(params, space, seed, sid) {
      // Set parameters unique to MD bead
      diameter_ = params->md_bead_diameter;
      mass_ = params->md_bead_mass;
    }
    ~MDBead() {}
    MDBead(const MDBead& that) : Bead(that) {}
    MDBead& operator=(MDBead const& that) {
      Bead::operator=(that); return *this;
    }
    virtual void Init();
    virtual void UpdatePosition();
    virtual void Integrate();
    virtual void UpdateKineticEnergy();
    virtual double const GetKineticEnergy();
};

class MDBeadSpecies : public Species<MDBead> {
  protected:
    void InitPotentials (system_parameters *params);
  public:
    MDBeadSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::md_bead);
      InitPotentials(params);
    }
    ~MDBeadSpecies() {}
    MDBeadSpecies(const MDBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }
};

#endif // _CYTOSCORE_MD_BEAD_H_
