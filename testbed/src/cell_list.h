// Most basic Cell List implementation I can think of
// Operates on base particles, with some given cutoff
#ifndef BUFFMD_CELLLIST_H_
#define BUFFMD_CELLLIST_H_

#include <vector>
#include <cmath>

#include "particle.h"

// Each actual cell has information about the ids of the particles
struct _cell {
    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int)*idxlist_.capacity();
        return mysize + vecsize;
    }
    
    int cell_id_; // ID of this cell
    int nparticles_; // Number of particles in this cell
    std::vector<int> idxlist_; // ID list of particles in this cell
};
typedef struct _cell cell_t;

class CellList {
public:
    
    CellList() {}
    ~CellList() {}
    
    void CreateCellList(int pN, double pRcut, double pSkin, double pBox[3]);
    void UpdateCellList(std::vector<particle*>* particles);
    void CheckCellList();
    
    unsigned long GetMemoryFootprint();
    
    // Simple getters
    int ncells() { return ncells_; }
    int npairs() { return npairs_; }
    int plist(int pairidx) {
        return plist_[pairidx];
    }
    
    // Operators
    cell_t* operator [](int cidx) {
        return &clist_[cidx];
    }
    

protected:
    
    // Inputs
    int nparticles_;
    double rcut_;
    double skin_;
    double rbuff_; // Buffer of rcut + skin
    double box_[3];
    
    // Computed quantities
    int ncells_;
    int npairs_;
    int nidx_;
    int T_[3];
    double S_[3];
    double boxby2_[3];
    
    // Vector of the cells themselves, and pair list
    std::vector<cell_t> clist_;
    std::vector<int> plist_;
    
    const int cellrat_ = 2;
    
};

#endif
