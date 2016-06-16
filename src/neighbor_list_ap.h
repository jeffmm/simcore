// Microcell list!

#ifndef _SIMCORE_NEIGHBOR_LIST_AP_H_
#define _SIMCORE_NEIGHBOR_LIST_AP_H_

#include "force_substructure_base.h"

struct _neighbor {
    int idx_;
    double value_;
};
typedef struct _neighbor neighbor_t;

typedef std::vector<neighbor_t> nl_list;

class NeighborListAP : public ForceSubstructureBase {
  public:

    NeighborListAP() {}
    virtual ~NeighborListAP() {
        delete[] neighbors_;
        printf("********\nNeighborList All Pairs stats\n");
        printf("\tNupdates: %d\n", n_updates_);
    }

    // Init should be the same
    // Load Simples should be the same
    virtual void CreateSubstructure(double pRcut);

    // Main update function
    void CheckNeighborList(bool pForceUpdate = false); // Checks for update if needed
    void UpdateNeighborList(); //Actually does the update
    void MinimumDistance(Simple *o1, Simple *o2, interactionmindist &idm);

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
