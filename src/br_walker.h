#ifndef _SIMCORE_BR_WALKER_H_
#define _SIMCORE_BR_WALKER_H_

#include "object.h"
#include "auxiliary.h"
#include "wca.h"

class BrWalker : public Simple {
  private:
    double diffusion_;
  public:
    BrWalker(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {
      diameter_=params->br_walker_diameter;
      SetDiffusion();
    }
    ~BrWalker() {}
    BrWalker(const BrWalker& that) : Simple(that) {}
    BrWalker& operator=(BrWalker const& that) {Simple::operator=(that); return *this;} 
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePosition();
    void UpdatePositionMP();
    void Init();
};

#include "species.h"
class BrWalkerSpecies : public Species<BrWalker> {
  protected:
    //void InitPotentials(system_parameters *params);

  public:
    BrWalkerSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_walker);
      //InitPotentials(params);
    }
    ~BrWalkerSpecies() {}
    BrWalkerSpecies(const BrWalkerSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _SIMCORE_BR_WALKER_H_

