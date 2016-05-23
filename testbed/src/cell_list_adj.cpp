// Cell list implementation

#include <cstring>
#include <cmath>
#include <iostream>
#include <cassert>

#include "cell_list_adj.h"

// First, generate the cells
// Independent of particle type, just need the best
// rcutoff to generate this
void
CellListAdj::CreateCellList(int pN, double pRcut, double pSkin, double pBox[3]){
    nparticles_ = pN;
    rcut_ = pRcut;
    memcpy(box_, pBox, 3*sizeof(double));
    
    SetRCut(rcut_, skin_);

    printf("Cell list has %dx%dx%d=%d cells of lengths {%.2f,%.2f,%.2f} "
           "with %d particles/cell.\n", T_[0], T_[1], T_[2], ncells_, S_[0], S_[1], S_[2],
           nidx_);
}

// Update Rcut, which requires rebuilding and checking lots of stuff
void
CellListAdj::SetRCut(double pRcut, double pSkin) {
    rcut_ = pRcut;
    skin_ = pSkin;
    
    // Recompute the number of cells on a side
    for (int i = 0; i < 3; ++i) {
        boxby2_[i] = 0.5*box_[i];
        T_[i] = floor(cellrat_ * box_[i] / (rcut_ + skin_));
        if (T_[i] < 3) {
            T_[i] = 3;
        }
    }
    
    // Compute S (side length) from this
    for (int i = 0; i < 3; ++i) {
        S_[i] = box_[i] / T_[i];
    }
    
    // Compute various needed quantities
    ncells_ = T_[0] * T_[1] * T_[2];
    // How many particles do we store in a cell?
    unsigned int density = 2*nparticles_ / ncells_ + 2;
    density = ((density/2)+1) * 2;
    nidx_ = (int)buffmd::nextpow2(density);
    if(nidx_ < 64) {
        nidx_ = 64;
    }
    
    clist_.clear();
    clist_.resize(ncells_);
    
    // Allocate index list within cells
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].cell_id_ = i;
        clist_[i].idxlist_.clear();
        clist_[i].idxlist_.resize(nidx_);
    }
    
    // Add the allowed adjacent cells to this one
    for (int cz = 0; cz < T_[2]; ++cz) {
        for (int cy = 0; cy < T_[1]; ++cy) {
            for (int cx = 0; cx < T_[0]; ++cx) {
                // Get the linear index of this cell
                int cidx = buffmd::cell_vec_to_linear(cx, cy, cz, T_);
                int adj_cell_id = 0;
                
                for (int dz = -1; dz <= 1; ++dz) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        for (int dx = -1; dx <= 1; ++dx) {
                            int rx = cx + dx;
                            int ry = cy + dy;
                            int rz = cz + dz;
                            
                            while(rx < 0) {
                                rx += T_[0];
                            }
                            rx = rx % T_[0];
                            
                            while(ry < 0) {
                                ry += T_[1];
                            }
                            ry = ry % T_[1];
                            
                            while(rz < 0) {
                                rz += T_[2];
                            }
                            rz = rz % T_[2];
                            
                            int cjdx = buffmd::cell_vec_to_linear(rx, ry, rz, T_);
                            clist_[cidx].adj_cell_ids_[adj_cell_id] = cjdx;
                            adj_cell_id++;
                        }
                    }
                }
                
            }
        }
    }
}


// Update the cell list
void
CellListAdj::UpdateCellList(std::vector<particle*>* particles) {
    int midx = 0;
    int cx, cy, cz;
    
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        clist_[cidx].nparticles_ = 0;
    }
    
    for (int i = 0; i < nparticles_; ++i) {
        int idx;
        
        auto p = (*particles)[i];
        cx = floor((buffmd::pbc(p->x[0], boxby2_[0], box_[0]) + boxby2_[0])/S_[0]);
        cy = floor((buffmd::pbc(p->x[1], boxby2_[1], box_[1]) + boxby2_[1])/S_[1]);
        cz = floor((buffmd::pbc(p->x[2], boxby2_[2], box_[2]) + boxby2_[2])/S_[2]);
        int cidx = buffmd::cell_vec_to_linear(cx,cy,cz, T_);
        //printf("p{%d}(%f,%f,%f) -> (%d,%d,%d) -> cell{%d}\n",
        //       i, p->x[0], p->x[1], p->x[2], cx, cy, cz, cidx);
        
        idx = clist_[cidx].nparticles_;
        clist_[cidx].idxlist_[idx] = i;
        ++idx;
        clist_[cidx].nparticles_ = idx;
        if (idx > midx) midx = idx;
    }
    
    // Check for overflow, shouldn't happen by far!
    if (midx > nidx_) {
        printf("Overflow in cell list: %d/%d particles/cells\n", midx, nidx_);
        exit(1);
    }
}


// Get the memory used by this class
unsigned long
CellListAdj::GetMemoryFootprint() {
    auto mysize = sizeof(*this);
    auto vecsize = clist_.capacity()*sizeof(CellAdj*);
    return mysize + vecsize;
}


// Check the cell list for consistency
void
CellListAdj::CheckCellList() {
    // Check the consistency of this new cell list
    CellAdj* memorycell = &clist_[0];
    printf("\t{CellListAdj size: %.1fkb}, {Cell size: %.1fb}, "
           "{total size: %.1fkb}\n",
           (float)(GetMemoryFootprint())/1024,
           (float)memorycell->GetMemoryFootprint(),
           (float)(GetMemoryFootprint() + ncells_*memorycell->GetMemoryFootprint())/1024);
    int runningtot = 0;
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        CellAdj* cell1 = &clist_[cidx];
        assert(cell1->cell_id_ == cidx);
        int cx[3];
        buffmd::cell_linear_to_vec(cidx, T_, cx);
        //printf("Cell(%d){%d,%d,%d} -> [%d]\n", cidx, cx[0], cx[1], cx[2], cell1->nparticles_);
        runningtot += cell1->nparticles_;
    }
    printf("Exact middle cell adjacent members: \n");
    int midcidx = buffmd::cell_vec_to_linear(0,0,0,T_);
    printf("Cell %d adjacent members = \n\t{", midcidx);
    for (int i = 0; i < 27; ++i) {
        printf("%d,", clist_[midcidx].adj_cell_ids_[i]);
    }


    //printf("Pair list: \n");
    //for (int pairidx = 0; pairidx < npairs_; ++pairidx) {
    //    cell_t* c1 = &clist_[2*pairidx];
    //    cell_t* c2 = &clist_[2*pairidx+1];
    //    printf("\t [%d] -> [%d]\n", c1->cell_id_, c2->cell_id_);
    //}
    assert(runningtot == nparticles_);
}


