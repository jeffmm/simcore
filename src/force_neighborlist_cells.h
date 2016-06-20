// Neighbor List cells force calculation

#ifndef _SIMCORE_FORCE_NEIGHBORLIST_CELLS_H_
#define _SIMCORE_FORCE_NEIGHBORLIST_CELLS_H_

#include "force_base.h"
#include "neighbor_list_cells.h"

class ForceNeighborListCells : public ForceBase {
  public:

    ForceNeighborListCells() {
        name_ = "ForceNeighborListCells";
    }
    virtual ~ForceNeighborListCells() {}

    // Override these functions, need special stuff
    virtual void Init(space_struct* pSpace, double pSkin);
    virtual void printSpecifics();
    virtual void dump();
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    NeighborListCells neighbor_list_;

};

#endif
