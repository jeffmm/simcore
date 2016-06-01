#ifndef _CYTOSCORE_MD_BEAD_H_
#define _CYTOSCORE_MD_BEAD_H_

#include "species.h"
#include "bead.h"
#include "auxiliary.h"
#include "test_potential.h"
#include "lennard_jones_12_6.h"

class MDBead : public Bead {
  protected:
    double prev_force_[3];
    double velocity_[3];
    double mass_;
  public:
    MDBead(system_parameters *params, space_struct *space, long seed, SID sid) : Bead(params, space, seed, sid) {
      diameter_ = params->md_bead_diameter;
      mass_ = params->md_bead_mass;
    }
    ~MDBead() {}
    MDBead(const MDBead& that) : Bead(that) {}
    MDBead& operator=(MDBead const& that) {Bead::operator=(that); return *this;} 
    void Init();
    void UpdatePosition();
    void Integrate();
    void UpdateKineticEnergy();
    double const GetKineticEnergy();
};

class MDBeadSpecies : public Species<MDBead> {
  protected:
    void InitPotentials () {
      AddPotential(SID::md_bead, SID::md_bead, new LJ126(1.0,1.0,space_, 10));
      AddPotential(SID::md_bead, SID::brownian_dimer, new TestPotential(space_, 10));
      AddPotential(SID::md_bead, SID::brownian_bead, new TestPotential(space_, 10));
    }

  public:
    MDBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::md_bead);
      InitPotentials();
    }
    ~MDBeadSpecies() {}
    MDBeadSpecies(const MDBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    //double GetTotalEnergy() {
      //double en=0;
      //for (auto it=members_.begin(); it!= members_.end(); ++it)
        //en += (*it)->GetEnergy();
      //return en;
    //}

};

#endif // _CYTOSCORE_MD_BEAD_H_
