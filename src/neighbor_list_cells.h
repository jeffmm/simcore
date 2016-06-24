// Microcell list!

#ifndef _SIMCORE_NEIGHBOR_LIST_CELLS_H_
#define _SIMCORE_NEIGHBOR_LIST_CELLS_H_

#include "force_substructure_base.h"
#include "neighbor_list_generic.h"
#include "adj_cell_list.h"

class NeighborListCells : public ForceSubstructureBase {
  public:

    NeighborListCells() {
        name_ = "NeighborListCells";
    }
    virtual ~NeighborListCells() {
        delete[] neighbors_;
        printf("********\nNeighborList Cells stats\n");
        printf("\tNupdates: %d, Skin:%2.2f\n", n_updates_, skin_);
    }

    // Override init for the cell list
    virtual void Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin);
    virtual void LoadFlatSimples(std::vector<Simple*> pSimples);
    virtual void print();
    virtual void dump();
    virtual void CreateSubstructure(double pRcut);

    // Main update function
    void CheckNeighborList(bool pForceUpdate = false); // Checks for update if needed
    void UpdateNeighborList(); //Actually does the update

    // Local cells update function
    void CellsUpdate();

    // Getters
    int GetNUpdates() { return n_updates_; }
    const nl_list* GetNeighbors() { return neighbors_; }

  protected:

    // Cell list quantities (computed)
    double rcut2_;
    double skin2_;
    double rcs2_;
    double half_skin2_;
    double boxby2_[3];

    bool nl_update_;
    int n_updates_;

    // The neighbor list itself
    nl_list* neighbors_;

    // Our underlying cell list
    AdjCellList cell_list_;

};

#endif
