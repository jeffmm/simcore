// Most basic Cell List implementation I can think of
// Operates on base particles, with some given cutoff
#ifndef BUFFMD_CELLLIST_H_
#define BUFFMD_CELLLIST_H_

#include <vector>
#include <cmath>

#include "particle.h"

// Each actual cell has information about the ids of the particles
class CellAdj {
public:
    
    CellAdj() {}
    ~CellAdj() {
        //delete(idxlist_);
    }

    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int)*idxlist_.capacity();
        return mysize + vecsize;
    }

    int cell_id_; // cell id
    int nparticles_; // number of particles in cell
    int adj_cell_ids_[13]; // Adjacent cell ids (including this one)
    std::vector<int> idxlist_;
    //int* idxlist_; // ID list of particles in cell
};

class CellListAdj {
public:
    
    CellListAdj() {}
    ~CellListAdj() {}
    
    void CreateCellList(int pN, double pRcut, double pBox[3]);
    void UpdateCellList(std::vector<particle*>* particles);
    void CheckCellList();
    
    unsigned long GetMemoryFootprint();
    
    // Simple getters
    int ncells() { return ncells_; }
    
    // Operators
    CellAdj* operator [](int cidx) {
        return &clist_[cidx];
    }
    

protected:
    
    // Inputs
    int nparticles_;
    double rcut_;
    double box_[3];
    
    // Computed quantities
    int ncells_;
    int npairs_;
    int nidx_;
    int T_[3];
    double S_[3];
    double boxby2_[3];
    
    // Vector of the cells themselves, and pair list
    std::vector<CellAdj> clist_;
    
    // Consts
    const int cellrat_ = 2;
    const int allowed_adj_[13][3] = {{1,0,0},{-1,-1,0},{0,-1,0},{1,-1,0},{-1,1,-1},{0,1,-1},{1,1,-1},
                                     {-1,0,-1},{0,0,-1},{1,0,-1},{-1,-1,-1},{0,-1,-1},{1,-1,-1}};
    
};

#endif
