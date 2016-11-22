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
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePosition();
    void UpdatePositionMP();
};

#include "species.h"
class BrBeadSpecies : public Species<BrBead> {
  protected:
    //void InitPotentials(system_parameters *params);

  public:
    BrBeadSpecies() : Species() {
      SetSID(SID::br_bead);
    }
    BrBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_bead);
    }
    void Init() {
      Species::Init();
    }

    void Configurator();
};

#endif // _SIMCORE_BR_BEAD_H_

