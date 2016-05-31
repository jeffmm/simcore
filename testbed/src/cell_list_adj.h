// Most basic Cell List implementation I can think of
// Operates on base particles, with some given cutoff
#ifndef BUFFMD_CELLLISTADJ_H_
#define BUFFMD_CELLLISTADJ_H_

#include <vector>
#include <cmath>

#include "particle.h"

// Each actual cell has information about the ids of the particles
class CellAdj {
public:
    
    CellAdj() {}
    ~CellAdj() {}

    unsigned long GetMemoryFootprint() {
        auto mysize = sizeof(*this);
        auto vecsize = sizeof(int)*idxlist_.capacity();
        return mysize + vecsize;
    }

    int cell_id_; // cell id
    int nparticles_; // number of particles in cell
    int adj_cell_ids_[27]; // Adjacent cell ids (including this one)
    std::vector<int> idxlist_;
};

class CellListAdj {
public:
    
    CellListAdj() {}
    ~CellListAdj() {}
    
    void CreateCellList(int pN, double pRcut, double pSkin, double pBox[3]);
    void UpdateCellList(std::vector<particle*>* particles);
    void CheckCellList();
    void dump();
    
    void SetRCut(double pRcut, double pSkin);
    
    unsigned long GetMemoryFootprint();
    
    // Simple getters
    int ncells() { return ncells_; }
    std::vector<int>* pidtocid() {
        return &p_c_;
    }
    
    // Operators (getters)
    CellAdj* operator [](int cidx) {
        return &clist_[cidx];
    }
    

protected:
    
    // Inputs
    int nparticles_;
    double rcut_;
    double skin_;
    double rbuff_;
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
    std::vector<int> p_c_;
    
};

#endif
