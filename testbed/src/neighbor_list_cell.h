// Main neighbor list implementation, should be derived from

// NOTE: Implements all pairs search by default!!!!

#ifndef BUFFMD_NEIGHBOR_LIST_CELL_H_
#define BUFFMD_NEIGHBOR_LIST_CELL_H_

#include <vector>

#include "particle.h"
#include "cell_list_adj.h"
#include "neighbor_list.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

class NeighborListCell : public NeighborList {
public:
    
    NeighborListCell() {}
    virtual ~NeighborListCell();
   
    virtual void print();
    virtual void dump();
    virtual void CreateNeighborList(int pN, double pRcut, double pSkin, double pBox[3]);
    virtual void UpdateNeighborList(std::vector<particle*>* particles);
    
    void CellUpdate(std::vector<particle*>* particles);
    void SetNThreads(int pNthreads) { nthreads_ = pNthreads; }
    
    
protected:
    
    // We inherit a ton of stuff from neighbor_list.h
    // But not cells
    CellListAdj cell_list_;
    
    int nthreads_;
};

#endif
