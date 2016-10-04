#ifndef _SIMCORE_POTENTIAL_MANAGER_V2_H_
#define _SIMCORE_POTENTIAL_MANAGER_V2_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "potential_base.h"
#include "species.h"

#include "helpers.h"

class PotentialManagerV2 {

  protected:
    int npots_;

    space_struct *space_;
    al_set *anchors_;
    std::vector<SpeciesBase*> *species_;

    rfh::factory pot_factory_;

    std::vector<PotentialBase*> potentials_;

  public:
    PotentialManagerV2() {}
    ~PotentialManagerV2() {}
    
    void Init(std::vector<SpeciesBase*> *pSpecies, space_struct *pSpace, al_set *pAnchors);
    void RegisterPotentials();
    void Print();

    int AddPotential(YAML::Node *subnode);
    PotentialBase *GetPotential(int idx) {return potentials_[idx];}
};

#endif // _SIMCORE_POTENTIAL_MANAGER_V2_H_

