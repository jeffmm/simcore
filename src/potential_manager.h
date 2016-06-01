#ifndef _CYTOSCORE_POTENTIAL_MANAGER_H_
#define _CYTOSCORE_POTENTIAL_MANAGER_H_

#include "auxiliary.h"
#include "potential_base.h"
#include "test_potential.h"

class PotentialManager {

  protected:
    potential_map potentials_;
  public:
    PotentialManager() {}
    ~PotentialManager() {
      potentials_.clear();
    }
    //void AddInteractions(std::vector<Interaction> interactions) {
      //for (auto it=interactions.begin(); it!= interactions.end(); ++it)
        //AddPotential(it->GetSID1(), it->GetSID2(), it->GetPotential());
    //}
    void AddPotential(SID sid1, SID sid2, PotentialBase *pot) {
      sid_pair key1 = std::make_pair(sid1, sid2);
      if (potentials_.count(key1)) return;
      sid_pair key2 = std::make_pair(sid2, sid1);
      if (potentials_.count(key2)) return;
      potentials_[key1] = pot;
    }
    PotentialBase * GetPotential(SID sid1, SID sid2) {
      sid_pair key1 = std::make_pair(sid1, sid2);
      if (potentials_.count(key1)) return potentials_[key1];
      sid_pair key2 = std::make_pair(sid2,sid1);
      if (potentials_.count(key2)) return potentials_[key2];
      return NULL;
    }
    void Print() {
      for (auto pot=potentials_.begin(); pot!=potentials_.end(); ++pot)
      {
        printf("{%d,%d} : ",(int) pot->first.first, (int) pot->first.second);
        pot->second->Print();
      }
    }
};

#endif // _CYTOSCORE_POTENTIAL_MANAGER_H_

