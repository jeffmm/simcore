#ifndef _CYTOSCORE_BROWNIAN_BEAD_H_
#define _CYTOSCORE_BROWNIAN_BEAD_H_

#include "bead.h"
#include "auxiliary.h"

class BrownianBead : public Bead {
  private:
    //void SetDiffusion();
    double diffusion_;
  public:
    BrownianBead(system_parameters *params, space_struct *space, long seed, unsigned int const sid) : Bead(params, space, seed, sid) {
      diameter_=params->br_bead_diameter;
      SetDiffusion();
    }
    ~BrownianBead() {}
    BrownianBead(const BrownianBead& that) : Bead(that) {}
    BrownianBead& operator=(BrownianBead const& that) {Bead::operator=(that); return *this;} 
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePosition();
};

#include "species.h"
class BrownianBeadSpecies : public Species<BrownianBead> {

  public:
    BrownianBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {}
    ~BrownianBeadSpecies() {}
    BrownianBeadSpecies(const BrownianBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _CYTOSCORE_BROWNIAN_BEAD_H_

