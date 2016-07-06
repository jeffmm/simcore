// Microcell list!

#ifndef _SIMCORE_NEIGHBOR_LIST_AP_H_
#define _SIMCORE_NEIGHBOR_LIST_AP_H_

#include "force_substructure_base.h"
#include "neighbor_list_generic.h"

class NeighborListAP : public ForceSubstructureBase {
  public:

    NeighborListAP() {
        name_ = "NeighborListAP";
    }
    virtual ~NeighborListAP() {
        delete[] neighbors_;
        printf("********\nNeighborList All Pairs stats\n");
        printf("\tNupdates: %d, Skin:%2.2f\n", n_updates_, skin_);
    }

    // Init should be the same
    // Load Simples should be the same
    virtual void print();
    virtual void dump();
    virtual void CreateSubstructure(double pRcut);

    // Main update function
    void CheckNeighborList(bool pForceUpdate = false); // Checks for update if needed
    void UpdateNeighborList(); //Actually does the update

    // Local all pairs update function
    void AllPairsUpdate();

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

};

#endif
