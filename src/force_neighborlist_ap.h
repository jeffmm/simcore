// Neighbor List all pairs force calculation

#ifndef _SIMCORE_FORCE_NEIGHBORLIST_AP_H_
#define _SIMCORE_FORCE_NEIGHBORLIST_AP_H_

#include "force_base.h"
#include "neighbor_list_ap.h"

class ForceNeighborListAP : public ForceBase {
  public:

    ForceNeighborListAP() {
        name_ = "ForceNeighborListAP";
    }
    virtual ~ForceNeighborListAP() {}

    // Override these functions, need special stuff
    virtual void Init(space_struct* pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);

    virtual void printSpecifics();
    virtual void dump();

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    NeighborListAP neighbor_list_;

};

#endif
