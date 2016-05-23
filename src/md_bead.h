#ifndef _CYTOSCORE_MD_BEAD_H_
#define _CYTOSCORE_MD_BEAD_H_

#include "species.h"
#include "bead.h"
#include "auxiliary.h"

class MDBead : public Bead {
  protected:
    double prev_force_[3];
    double velocity_[3];
    double mass_;
    double energy_;
  public:
    MDBead(system_parameters *params, space_struct *space, long seed) : Bead(params, space, seed) {
      diameter_ = params->md_bead_diameter;
      mass_ = params->md_bead_mass;
    }
    ~MDBead() {}
    MDBead(const MDBead& that) : Bead(that) {}
    MDBead& operator=(MDBead const& that) {Bead::operator=(that); return *this;} 
    void Init();
    void UpdatePosition();
    void Integrate();
    void UpdateEnergy();
    double const GetEnergy();
};

class MDBeadSpecies : public Species<MDBead> {

  public:
    MDBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {}
    ~MDBeadSpecies() {}
    MDBeadSpecies(const MDBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }

    double GetTotalEnergy() {
      double en=0;
      for (auto it=members_.begin(); it!= members_.end(); ++it)
        en += (*it)->GetEnergy();
      return en;
    }
};

#endif // _CYTOSCORE_MD_BEAD_H_
