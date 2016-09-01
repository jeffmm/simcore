#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "auxiliary.h"
#include "potential_base.h"
#include "species.h"
#include "test_potential.h"

#include "helpers.h"

#include <unordered_map>

class PotentialManager {

  protected:
    space_struct *space_;
    std::string fname_;
    int npots_;
    potential_map potentials_;
    // XXX Possibly just store the location of the potential base
    std::unordered_map<std::pair<unsigned int, unsigned int>, int, hashh::pair_hash> tethers_;
    std::vector<PotentialBase*> internal_potentials_;
    std::vector<SpeciesBase*> *species_;

    rfh::factory pot_factory_;

  public:
    PotentialManager() {}
    ~PotentialManager() {
      potentials_.clear();
    }
    
    void Init(std::vector<SpeciesBase*> *pSpecies, space_struct *pSpace, char *pFname);
    void RegisterPotentials();
    void ParsePotentials();

    double GetMaxRCut();

    void AddPotential(SID sid1, SID sid2, PotentialBase *pot);
    PotentialBase * GetPotential(SID sid1, SID sid2);
    PotentialBase * GetPotentialTether(unsigned int oid1, unsigned int oid2);
    void Print();
};

#endif // _SIMCORE_POTENTIAL_MANAGER_H_

