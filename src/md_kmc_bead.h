#ifndef _SIMCORE_MD_KMC_BEAD_H_
#define _SIMCORE_MD_KMC_BEAD_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"
#include "lennard_jones_12_6.h"
#include "wca.h"

class MDKMCBead : public Simple {
  protected:
    double mass_;
  public:
    MDKMCBead(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {
      // Set parameters unique to MD bead
      diameter_ = params->md_bead_diameter;
      mass_ = params->md_bead_mass;
    }
    ~MDKMCBead() {}
    MDKMCBead(const MDKMCBead& that) : Simple(that) {}
    MDKMCBead& operator=(MDKMCBead const& that) {
      Simple::operator=(that); return *this;
    }
    virtual void Init();
    virtual void UpdatePosition();
    virtual void UpdatePositionMP();
    virtual void Integrate();
    virtual void UpdateKineticEnergy();
    virtual double const GetKineticEnergy();
};

class MDKMCBeadSpecies : public Species<MDKMCBead> {
  protected:
    void InitPotentials (system_parameters *params);
  public:
    MDKMCBeadSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::md_bead);
      InitPotentials(params);
    }
    ~MDKMCBeadSpecies() {}
    MDKMCBeadSpecies(const MDKMCBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }
};

#endif // _SIMCORE_MD_KMC_BEAD_H_
