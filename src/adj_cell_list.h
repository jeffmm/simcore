// Adjacent cell list

#ifndef _SIMCORE_ADJ_CELL_LIST_H_
#define _SIMCORE_ADJ_CELL_LIST_H_

#include "force_substructure_base.h"

class Adjcell {
  public:

    Adjcell() {}
    ~Adjcell() {}

    // Helper for getting memory footprint
    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int)*idxlist_.capacity();
        return mysize + vecsize;
    }

    int cell_id_;
    int nparticles_;
    int* adj_cell_ids_;
    std::vector<int> idxlist_;
};

class AdjCellList : public ForceSubstructureBase {
  public:

    AdjCellList() {
        name_ = "AdjCellList";
    }
    virtual ~AdjCellList() {
        for (int i = 0; i < ncells_; ++i) {
            delete[] clist_[i].adj_cell_ids_;
        }
    }

    // Init should be the same
    // Load Simples should be the same
    virtual void print();
    virtual void CreateSubstructure(double pRcut);

    // Local function for updating cell list
    void UpdateCellList();

    // Simple getters
    int ncells() { return ncells_; }
    int nadj() { return nadj_; }
    std::vector<int>* pidtocid() { return &p_c_; }

    // Operators
    Adjcell* operator [](int cidx) { return &clist_[cidx]; }

  protected:

    // Cell list quantities (computed)
    int ncells_;
    int nadj_;
    int nidx_;
    int T_[3] = {1, 1, 1};
    double S_[3] = {0.0, 0.0, 0.0};;
    double boxby2_[3];

    // Needed for cells and pairs themselves
    std::vector<Adjcell> clist_;
    std::vector<int> p_c_; // particle to cell indexer

};

#endif
