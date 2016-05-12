// Cell list implementation v2

#ifndef TESTBD_CELL_LIST_V2_H_
#define TESTBD_CELL_LIST_V2_H_

#include <vector>
#include <memory>

#include "particle.h"

class CellTraditional {
public:
    CellTraditional() {
        memcpy(adjacent_cells_, adjacent_cells_ + 27*sizeof(int), 0);
    }
    ~CellTraditional() {}
    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int) * idxlist_.capacity();
        return mysize + vecsize;
    }
    
    int cell_id_; // ID of this cell
    int nparticles_; // Number of particles in this cell
    int adjacent_cells_[27]; // Always have 27-1 adjacent cells
    std::vector<int> idxlist_; // List of particle IDs in this cell
};

class CellListTraditional {
public:

    CellListTraditional() {}
    ~CellListTraditional() {}
    
    void CreateCellList(int pN, double pRcut, double pBox[3]);
    void UpdateCellList(std::vector<particle*>* particles);
    void CheckCellList();
    
    unsigned long GetMemoryFootprint();
    
    // Simple getters
    int ncells() { return ncells_; }
    

private:

    // Inputs
    int nparticles_; // Total number of cell list particles
    double rcut_; // Radial cutoff distance
    double box_[3]; // Unit cell size

    // Computed quantities
    int ncells_; // Total number of cells
    int nidx_; // Number of particles allowed in each cell
    int T_[3]; // number of cells per dimension
    double S_[3]; // Cell length per dimension
    double boxby2_[3]; // Half of the box size
    
    // Vector of the cells themselves
    std::vector<CellTraditional> clist_;
};

#endif
