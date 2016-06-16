// Microcell force computation

#ifndef _SIMCORE_FORCE_NEIGHBORLISTAP_H_
#define _SIMCORE_FORCE_NEIGHBORLISTAP_H_

#include "force_base.h"
#include "neighbor_list_ap.h"

class ForceNeighborListAP : public ForceBase {
  public:

    ForceNeighborListAP() {}
    virtual ~ForceNeighborListAP() {}

    // Override these functions, need special stuff
    virtual void Init(space_struct* pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    NeighborListAP neighbor_list_;

};

#endif
