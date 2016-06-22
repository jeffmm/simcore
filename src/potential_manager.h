#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "auxiliary.h"
#include "potential_base.h"
#include "test_potential.h"

#include "helpers.h"

class PotentialManager {

  protected:
    space_struct *space_;
    std::string fname_;
    int npots_;
    potential_map potentials_;

    rfh::factory pot_factory_;

  public:
    PotentialManager() {}
    ~PotentialManager() {
      potentials_.clear();
    }
    
    void Init(space_struct *pSpace, char *pFname);
    void RegisterPotentials();
    void ParsePotentials();

    double GetMaxRCut();

    void AddPotential(SID sid1, SID sid2, PotentialBase *pot);
    PotentialBase * GetPotential(SID sid1, SID sid2);
    void Print();
};

#endif // _SIMCORE_POTENTIAL_MANAGER_H_

