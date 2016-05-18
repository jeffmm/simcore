// Cell list implementation

#include <cstring>
#include <cmath>
#include <iostream>
#include <cassert>

#include "cell_list.h"

// First, generate the cells
// Independent of particle type, just need the best
// rcutoff to generate this
void
CellList::CreateCellList(int pN, double pRcut, double pBox[3]){
    double boxoffs[3];
    
    nparticles_ = pN;
    rcut_ = pRcut;
    memcpy(box_, pBox, 3*sizeof(double));
    
    // Compute the number of cells on a side
    for (int i = 0; i < 3; ++i) {
        boxby2_[i] = 0.5*box_[i];
        T_[i] = floor(cellrat_ * box_[i] / rcut_);
        if (T_[i] < 3) {
            T_[i] = 3;
        }
    }
    
    // Compute S (side length) from this
    for (int i = 0; i < 3; ++i) {
        S_[i] = box_[i] / T_[i];
        boxoffs[i] = boxby2_[i] - 0.5*S_[i];
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
    plist_.clear();
    plist_.resize(2*ncells_*ncells_);
    
    // Allocate index list within cells
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].cell_id_ = i;
        clist_[i].idxlist_.clear();
        clist_[i].idxlist_.resize(nidx_);
    }
    

    // Build the cell pair list
    npairs_ = 0;
    for (int cidx = 0; cidx < ncells_ - 1; ++cidx) {
        double x1, x2, y1, y2, z1, z2, rx, ry, rz;
        
        int cx1[3];
        buffmd::cell_linear_to_vec(cidx, T_, cx1);
        x1 = cx1[0]*S_[0] - boxoffs[0];
        y1 = cx1[1]*S_[1] - boxoffs[1];
        z1 = cx1[2]*S_[2] - boxoffs[2];
        
        for (int cjdx = cidx+1; cjdx < ncells_; ++cjdx) {
            int cx2[3];
            buffmd::cell_linear_to_vec(cjdx, T_, cx2);
            x2 = cx2[0]*S_[0] - boxoffs[0];
            y2 = cx2[1]*S_[1] - boxoffs[1];
            z2 = cx2[2]*S_[2] - boxoffs[2];
            
            rx = buffmd::pbc(x1 - x2, boxby2_[0], box_[0]);
            ry = buffmd::pbc(y1 - y2, boxby2_[1], box_[1]);
            rz = buffmd::pbc(z1 - z2, boxby2_[2], box_[2]);
            
            // Check the cells on a line that are too far apart
            if (fabs(rx) > rcut_ + S_[0]) continue;
            if (fabs(ry) > rcut_ + S_[1]) continue;
            if (fabs(rz) > rcut_ + S_[2]) continue;
            
            // Check for cells in a plane that are too far apart
            if (sqrt(rx*rx + ry*ry) > (rcut_ + sqrt(S_[0]*S_[0] + S_[1]*S_[1]))) continue;
            if (sqrt(rx*rx + rz*rz) > (rcut_ + sqrt(S_[0]*S_[0] + S_[2]*S_[2]))) continue;
            if (sqrt(ry*ry + rz*rz) > (rcut_ + sqrt(S_[1]*S_[1] + S_[2]*S_[2]))) continue;
            
            // Other cells too far apart in 3d
            if (sqrt(rx*rx + ry*ry + rz*rz) > (sqrt(S_[0]*S_[0] + S_[1]*S_[1] + S_[2]*S_[2]) + rcut_)) continue;
            
            // Cells are close enough, add
            plist_[2*npairs_    ] = cidx;
            plist_[2*npairs_ + 1] = cjdx;
            ++npairs_;
        }
    }
    printf("********\n");
    printf("Cell list has %dx%dx%d=%d cells of lengths {%.2f,%.2f,%.2f} "
           "with %d/%d pairs and %d particles/cell.\n", T_[0], T_[1], T_[2], ncells_, S_[0], S_[1], S_[2],
           npairs_, ncells_*(ncells_-1)/2, nidx_);
}


// Update the cell list
void
CellList::UpdateCellList(std::vector<particle*>* particles) {
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
CellList::GetMemoryFootprint() {
    auto mysize = sizeof(*this);
    auto vecsize = clist_.capacity()*sizeof(cell_t*);
    auto pvecsize = plist_.capacity()*sizeof(int*);
    return mysize + vecsize + pvecsize;
}


// Check the cell list for consistency
void
CellList::CheckCellList() {
    // Check the consistency of this new cell list
    cell_t* memorycell = &clist_[0];
    printf("\t{CellList size: %.1fkb}, {Cell size: %.1fb}, "
           "{total size: %.1fkb}\n",
           (float)(GetMemoryFootprint())/1024,
           (float)memorycell->GetMemoryFootprint(),
           (float)(GetMemoryFootprint() + ncells_*memorycell->GetMemoryFootprint())/1024);
    int runningtot = 0;
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        cell_t* cell1 = &clist_[cidx];
        assert(cell1->cell_id_ == cidx);
        int cx[3];
        buffmd::cell_linear_to_vec(cidx, T_, cx);
        //printf("Cell(%d){%d,%d,%d} -> [%d]\n", cidx, cx[0], cx[1], cx[2], cell1->nparticles_);
        runningtot += cell1->nparticles_;
    }
    //printf("Pair list: \n");
    //for (int pairidx = 0; pairidx < npairs_; ++pairidx) {
    //    cell_t* c1 = &clist_[2*pairidx];
    //    cell_t* c2 = &clist_[2*pairidx+1];
    //    printf("\t [%d] -> [%d]\n", c1->cell_id_, c2->cell_id_);
    //}
    assert(runningtot == nparticles_);
}


