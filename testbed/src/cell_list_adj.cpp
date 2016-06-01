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
    skin_ = pSkin;
    memcpy(box_, pBox, 3*sizeof(double));
    p_c_.clear();
    p_c_.resize(nparticles_);

    #if defined(_OPENMP)
    #pragma omp parallel
    {
        if(0 == omp_get_thread_num()) {
            nthreads_ = omp_get_num_threads();
        }
    }
    #else
    nthreads_ = 1;
    #endif
    
    SetRCut(rcut_, skin_);

    printf("Adj. Cell list has %dx%dx%d=%d cells of lengths {%.2f,%.2f,%.2f} "
           "with %d particles/cell (tot:%d).\n", T_[0], T_[1], T_[2], ncells_, S_[0], S_[1], S_[2],
           nidx_, nparticles_);
    printf("Done creating cell list adj.\n");
}

// Update Rcut, which requires rebuilding and checking lots of stuff
void
CellListAdj::SetRCut(double pRcut, double pSkin) {
    rcut_ = pRcut;
    skin_ = pSkin;
    rbuff_ = rcut_ + skin_;
    
    // Recompute the number of cells on a side
    for (int i = 0; i < 3; ++i) {
        boxby2_[i] = 0.5*box_[i];
        T_[i] = floor(box_[i] / rbuff_);
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
        
        idx = clist_[cidx].nparticles_;
        clist_[cidx].idxlist_[idx] = i;
        p_c_[i] = cidx; // Set the particle id->cell id
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


// Get the interaction pairs for this cell list?
void
CellListAdj::InteractionPairs(std::vector<std::pair<int, int>>* pPartPairs) {
    // Clear the vector, rebuild
    // XXX:
    // Right now this is awful, since we use an omp critical section, and we
    // really shouldn't need to, but because of vector memory allocs, we were
    // running into segfaults when trying to access it correctly!
    // DOESN'T WORK RIGHT NOW!!!
    //printf("Interaction pairs have serious problems, don't use them!\n");

    (*pPartPairs).clear();
    size_t *prefix;
   
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        int tid;
        #if defined(_OPENMP)
        tid = omp_get_thread_num();
        #else
        tid = 0;
        #endif

        #if defined(_OPENMP)
        #pragma omp single
        #endif
        {
            prefix = new size_t[nthreads_+1];
            prefix[0] = 0;
        }
        std::vector<std::pair<int,int>> vec_private;

        // Loop over all particles and build the list
        #if defined(_OPENMP)
        #pragma omp for schedule(runtime) nowait
        #endif
        for (int idx = 0; idx < nparticles_; ++idx) {
            // Get our cell
            int cidx = p_c_[idx];
            auto cell1 = clist_[cidx];
            // Loop over other cells (including us) in the block of cells
            for (int cjdx = 0; cjdx < 27; ++cjdx) {
                auto cell2 = clist_[cell1.adj_cell_ids_[cjdx]];
                // Loop over it's particles
                for (int jdx = 0; jdx < cell2.nparticles_; ++jdx) {
                    int jjdx = cell2.idxlist_[jdx];
                    // ONLY DO THE CALCULATION IF THE OTHER PID IS HIGHER!!!
                    if (jjdx > idx) {
                        vec_private.push_back(std::make_pair(idx, jjdx));
                    } // only do the calculation if the second id is higher
                } // cell2 particle list
            } // cells adjancent and equal to us
        } // pragma omp for schedule(runtime)
        prefix[tid+1] = vec_private.size();

        #if defined(_OPENMP)
        #pragma omp barrier
        #pragma omp single
        #endif
        {
            for (int i = 1; i < (nthreads_+1); ++i) prefix[i] += prefix[i-1];
            (*pPartPairs).resize((*pPartPairs).size() + prefix[nthreads_]);
        }
        std::copy(vec_private.begin(), vec_private.end(), (*pPartPairs).begin() + prefix[tid]);
    } // omp parallel

    delete[] prefix;
}


// Print out everything in gory detail
void
CellListAdj::dump() {
    printf("\n********\n");
    printf("Dumping Adj. Cell List!\n\n");
    printf("Adj. Cell list has %dx%dx%d=%d cells of lengths {%.2f,%.2f,%.2f} "
           "with %d particles/cell (tot:%d).\n", T_[0], T_[1], T_[2], ncells_, S_[0], S_[1], S_[2],
           nidx_, nparticles_);
    // Loop over the cells, and print their information
    printf("Printing cells!\n");
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        assert(cidx == clist_[cidx].cell_id_);
        auto cell1 = clist_[cidx];
        printf("Cell(%d) -> {n:%d}\n", cell1.cell_id_, cell1.nparticles_);
        // Print out the adjacent cell ids
        printf("\t adj: ");
        for (int i = 0; i < 27; ++i) {
            printf("%d,", cell1.adj_cell_ids_[i]);
        }
        printf("\n");
        printf("\t ids: ");
        for (int idx = 0; idx < cell1.nparticles_; ++idx) {
            printf("%d,", cell1.idxlist_[idx]);
        }
        printf("\n");
    }
    printf("\n********\n");
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
    printf("Exact bottom cell adjacent members: \n");
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


