// Cell force computation

#ifndef _SIMCORE_FORCE_CELL_H_
#define _SIMCORE_FORCE_CELL_H_

#include "force_base.h"
#include "adj_cell_list.h"

class ForceCell : public ForceBase {
  public:

    ForceCell() {}
    virtual ~ForceCell() {}

    // Override these functions, need special stuff
    virtual void Init(space_struct* pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    AdjCellList cell_list_;
};

#endif
