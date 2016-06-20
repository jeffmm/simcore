// Microcell list!

#ifndef _SIMCORE_MICROCELL_LIST_H_
#define _SIMCORE_MICROCELL_LIST_H_

#include "force_substructure_base.h"

class Microcell {
  public:

    Microcell() {}
    ~Microcell() {}

    // Helper for getting memory footprint
    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int)*idxlist_.capacity();
        return mysize + vecsize;
    }

    int cell_id_;
    int nparticles_;
    std::vector<int> idxlist_;
};

class MicrocellList : public ForceSubstructureBase {
  public:

    MicrocellList() {
        name_ = "MicrocellList";
    }
    virtual ~MicrocellList() {}

    // Init should be the same
    // Load Simples should be the same
    virtual void print();
    virtual void dump();
    virtual void CreateSubstructure(double pRcut);

    // Local function for updating cell list
    void UpdateCellList();

    // Simple getters
    int ncells() { return ncells_; }
    int npairs() { return npairs_; }
    int plist(int pairidx) { return plist_[pairidx]; }
    std::vector<int>* pidtocid() { return &p_c_; }

    // Operators
    Microcell* operator [](int cidx) { return &clist_[cidx]; }

  protected:

    // Cell list quantities (computed)
    int ncells_;
    int npairs_;
    int nidx_;
    int T_[3] = {1, 1, 1};
    double S_[3] = {0.0, 0.0, 0.0};;
    double boxby2_[3];

    // Needed for cells and pairs themselves
    std::vector<Microcell> clist_;
    std::vector<int> plist_;
    std::vector<int> p_c_; // particle to cell indexer

    const int cellrat_ = 2.0;
};

#endif
