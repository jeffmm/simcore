#ifndef _SIMCORE_BR_BEAD_H_
#define _SIMCORE_BR_BEAD_H_

#include "object.h"
#include "auxiliary.h"
#include "wca.h"

class BrBead : public Simple {
  private:
    double diffusion_;
  public:
    BrBead(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {
      diameter_=params->br_bead_diameter;
      SetDiffusion();
    }
    ~BrBead() {}
    BrBead(const BrBead& that) : Simple(that) {}
    BrBead& operator=(BrBead const& that) {Simple::operator=(that); return *this;} 
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePosition();
};

#include "species.h"
class BrBeadSpecies : public Species<BrBead> {
  protected:
    void InitPotentials(system_parameters *params);

  public:
    BrBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_bead);
      InitPotentials(params);
    }
    ~BrBeadSpecies() {}
    BrBeadSpecies(const BrBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _SIMCORE_BR_BEAD_H_

