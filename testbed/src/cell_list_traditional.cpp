// Implementation for traditional cell lists

#include <cstring>
#include <cmath>
#include <iostream>
#include <cassert>

#include "cell_list.h"
#include "helpers.h"

void
CellListTraditional::CreateCellList(int pN, double pRcut, double pBox[3]){
    nparticles_ = pN;
    rcut_ = pRcut;
    memcpy(box_, pBox, 3*sizeof(double));

    // Compute the number of cells on a side
    for (int i = 0; i < 3; ++i) {
        boxby2_[i] = box_[i];
        T_[i] = ceil(box_[i] / rcut_);
    }
    // Minimum of 3 cells on a side
    for (int i = 0; i < 3; ++i) {
        if (T_[i] < 3) {
            T_[i] = 3;
        }
    }
    // Compute S from this (in case under the limit of 3 per side)
    for (int i = 0; i < 3; ++i) {
        S_[i] = box_[i] / T_[i];
    }
    
    // Compute various needed quantities
    ncells_ = T_[0] * T_[1] * T_[2];
    // How many particles per cell should we allow?  Upper bound of 2^density or 64
    unsigned int density = 2*nparticles_ / ncells_ + 2;
    density = ((density/2)+1) * 2;
    nidx_ = (int)buffmd::nextpow2(density);
    if(nidx_ < 64) {
        nidx_ = 64;
    }
    
    clist_.clear();
    clist_.resize(ncells_);
    
    // Allocate the index list within cells
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].idxlist_.clear();
        clist_[i].idxlist_.resize(nidx_);
    }
    
    // For each cell, find it's neighboring cells (plus periodic boundary conditions!)
    for (int cz = 0; cz < T_[2]; ++cz) {
        for (int cy = 0; cy < T_[1]; ++cy) {
            for (int cx = 0; cx < T_[0]; ++cx) {
                // Calculate the scalar index of the cell
                int cidx = buffmd::cell_vec_to_linear((int []){cx, cy, cz}, T_);
                CellTraditional* cell1 = &clist_[cidx];
                // Who are all of my neighbors?
                int i_neighb = 0;
                for (int z = -1; z <= 1; ++z) {
                    for (int y = -1; y <= 1; ++y) {
                        for (int x = -1; x <= 1; ++x) {
                            int rx = cx + x;
                            if (rx < 0)
                                rx += T_[0];
                            rx = rx % T_[0];
                            
                            int ry = cy + y;
                            if (ry < 0)
                                ry += T_[1];
                            ry = ry % T_[1];
                            
                            int rz = cz + z;
                            if (rz < 0)
                                rz += T_[2];
                            rz = rz % T_[2];

                            int cidx2 = buffmd::cell_vec_to_linear((int[]){rx,ry,rz}, T_);
                            cell1->adjacent_cells_[i_neighb++] = cidx2;
                        }
                    }
                }
                // Do not sort the resulting adjacent array for cache purposes
            }
        }
    }
}

void
CellListTraditional::UpdateCellList(std::vector<std::shared_ptr<particleBase>>& particles) {
    // Reset the cell list
    int cx, cy, cz;

    // Reset number of particles in each box
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].nparticles_ = 0;
    }
    
    int midx = 0;
    for (int i = 0; i < nparticles_; ++i) {
        int idx;
        
        auto p = particles[i];
        cx = floor((buffmd::pbc(p->x[0], boxby2_[0], box_[0]) + boxby2_[0])/S_[0]);
        cy = floor((buffmd::pbc(p->x[1], boxby2_[1], box_[1]) + boxby2_[1])/S_[1]);
        cz = floor((buffmd::pbc(p->x[2], boxby2_[2], box_[2]) + boxby2_[2])/S_[2]);
        int cidx = buffmd::cell_vec_to_linear((int[]){cx,cy,cz}, T_);
        printf("p{%d}(%f,%f,%f) -> (%d,%d,%d) -> cell{%d}\n",
               i, p->x[0], p->x[1], p->x[2], cx, cy, cz, cidx);
    }
}
