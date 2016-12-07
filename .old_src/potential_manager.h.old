#ifndef _SIMCORE_POTENTIAL_MANAGER_H_
#define _SIMCORE_POTENTIAL_MANAGER_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "potential_base.h"
#include "species.h"
#include "test_potential.h"

#include "helpers.h"

#include <unordered_map>

class PotentialManager {

  protected:
    space_struct *space_;
    al_set *anchors_;
    std::string fname_;
    int npots_;
    potential_map potentials_;
    std::vector<PotentialBase*> potential_vec_;
    std::vector<std::string> potential_vec_names_types_;
    std::vector<PotentialBase*> potential_vec_tethers_;
    std::vector<std::string> potential_vec_names_tethers_;
    std::unordered_map<std::pair<unsigned int, unsigned int>, PotentialBase*, hashh::pair_hash> internal_potentials_;
    std::unordered_map<std::pair<unsigned int, unsigned int>, PotentialBase*, hashh::pair_hash> tethers_;
    std::map<SID, PotentialBase*> boundaries_;
    std::vector<SpeciesBase*> *species_;
    std::vector<SID> internal_sids_;

    rfh::factory pot_factory_;

  public:
    PotentialManager() {}
    ~PotentialManager() {
      potentials_.clear();
    }
    
    void Init(std::vector<SpeciesBase*> *pSpecies, space_struct *pSpace, al_set *pAnchors, char *pFname);
    void RegisterPotentials();
    void ParsePotentials();

    double GetMaxRCut();

    void AddPotentialExternal(SID sid1, SID sid2, PotentialBase *pot);
    void AddPotentialInternal(unsigned int oid1, unsigned int oid2, PotentialBase *pot);
    void AddPotentialTether(SID sid1, SID sid2, PotentialBase *pot);
    PotentialBase * GetPotentialExternal(SID sid1, SID sid2);
    PotentialBase * GetPotentialInternal(unsigned int oid1, unsigned int oid2);
    PotentialBase * GetPotentialTether(unsigned int oid1, unsigned int oid2);
    PotentialBase * GetPotentialBoundary(SID sid);
    std::vector<PotentialBase*> GetAllPotentials();
    std::vector<SID> *GetInternalSIDs() {return &internal_sids_;}
    void Print();
};

#endif // _SIMCORE_POTENTIAL_MANAGER_H_

