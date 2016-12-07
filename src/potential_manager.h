#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "potential_base.h"
#include "species.h"

#include "helpers.h"

class PotentialManager {

  protected:
    int npots_;

    space_struct *space_;
    al_set *anchors_;

    rfh::factory pot_factory_;

    std::vector<PotentialBase*> potentials_;

  public:
    PotentialManager() {}
    ~PotentialManager() {}
    
    void Init(space_struct *pSpace, al_set *pAnchors);
    void RegisterPotentials();
    void Print();

    int AddPotential(YAML::Node *subnode);
    PotentialBase *GetPotential(int idx) {return potentials_[idx];}
};

#endif // _SIMCORE_POTENTIAL_MANAGER_H_

