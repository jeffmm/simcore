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
    //double velocity_[3];
    double mass_;
  public:
    MDBead(system_parameters *params, space_struct *space, long seed, SID sid) : Bead(params, space, seed, sid) {
      diameter_ = params->md_bead_diameter;
      mass_ = params->md_bead_mass;
    }
    ~MDBead() {}
    MDBead(const MDBead& that) : Bead(that) {}
    MDBead& operator=(MDBead const& that) {Bead::operator=(that); return *this;} 
    virtual void Init();
    virtual void UpdatePosition();
    virtual void Integrate();
    virtual void UpdateKineticEnergy();
    virtual double const GetKineticEnergy();
};

class MDBeadSpecies : public Species<MDBead> {
  protected:
    void InitPotentials (system_parameters *params) {
      AddPotential(SID::md_bead, SID::md_bead, new LJ126(params->lj_epsilon,params->md_bead_diameter,space_, 2.5*params->md_bead_diameter));
      AddPotential(SID::md_bead, SID::brownian_dimer, new TestPotential(space_, 10));
      AddPotential(SID::md_bead, SID::brownian_bead, new TestPotential(space_, 10));
    }

  public:
    MDBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
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
      //double temp[3] = {-4,0,0};
      //members_[0]->SetPosition(temp);
      //temp[0] = 4;
      //members_[1]->SetPosition(temp);
      //temp[0] = 1;
      //members_[0]->SetVelocity(temp);
      //temp[0] = -1;
      //members_[1]->SetVelocity(temp);
      //temp[0] = -4 - 0.0001;
      //members_[0]->SetPrevPosition(temp);
      //temp[0] = 4 + 0.0001;
      //members_[1]->SetPrevPosition(temp);
      //for (auto it=members_.begin(); it!=members_.end(); ++it) {
        //(*it)->UpdatePeriodic();
      //}
    //double GetTotalEnergy() {
      //double en=0;
      //for (auto it=members_.begin(); it!= members_.end(); ++it)
        //en += (*it)->GetEnergy();
      //return en;
    //}

};

#endif // _CYTOSCORE_MD_BEAD_H_
