// implementation for adjacent cell list

#include <chrono>

#include "adj_cell_list.h"

#include "minimum_distance.h"

void
AdjCellList::CreateSubstructure(double pRcut) {
    printf("AdjCellList::CreateSubstructure begin\n");
    auto start = std::chrono::steady_clock::now();

    rcut_ = pRcut;
    rbuff_ = rcut_ + skin_;
    p_c_.clear();
    p_c_.resize(nparticles_);

    // Compute the number of cells on a side
    for (int i = 0; i < ndim_; ++i) {
        boxby2_[i] = 0.5*box_[i];
        T_[i] = floor(box_[i] / rbuff_);
        if (T_[i] < 3) {
            T_[i] = 3;
        }
    }

    // Compute the side length
    for (int i = 0; i < ndim_; ++i) {
        S_[i] = box_[i] / T_[i];
    }

    // Compute various needed quantities
    ncells_ = 1;
    for (int i = 0; i < ndim_; ++i) {
        ncells_ *= T_[i];
    }
    // How many particles should we store in a cell?
    unsigned int density = 2*nparticles_ / ncells_ + 2;
    density = ((density/2)+1) * 2;
    nidx_ = (int)cytohelpers::nextpow2(density);
    if (nidx_ < 64) nidx_ = 64;

    clist_.clear();
    clist_.resize(ncells_);

    // Allocate the list index within cells
    nadj_ = pow(3, ndim_);
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].cell_id_ = i;
        clist_[i].adj_cell_ids_ = new int[nadj_];
        clist_[i].idxlist_.clear();
        clist_[i].idxlist_.resize(nidx_);
    }

    if (ndim_ == 3) {
        // Add the allowed adjacent cells to this one
        for (int cz = 0; cz < T_[2]; ++cz) {
            for (int cy = 0; cy < T_[1]; ++cy) {
                for (int cx = 0; cx < T_[0]; ++cx) {
                    // Get the linear index of this cell
                    int cidx = cytohelpers::cell_vec_to_linear(cx, cy, cz, T_);
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
                                
                                int cjdx = cytohelpers::cell_vec_to_linear(rx, ry, rz, T_);
                                clist_[cidx].adj_cell_ids_[adj_cell_id] = cjdx;
                                adj_cell_id++;
                            }
                        }
                    }
                    
                }
            }
        }
    } else if (ndim_ == 2) {
        for (int cy = 0; cy < T_[1]; ++cy) {
            for (int cx = 0; cx < T_[0]; ++cx) {
                // Get the linear index of the cell (cz = 0)
                int cidx = cytohelpers::cell_vec_to_linear(cx, cy, 0, T_);
                int adj_cell_id = 0;

                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int rx = cx + dx;
                        int ry = cy + dy;

                        while (rx < 0) {
                            rx += T_[0];
                        }
                        rx = rx % T_[0];

                        while (ry < 0) {
                            ry += T_[1];
                        }
                        ry = ry % T_[1];

                        int cjdx = cytohelpers::cell_vec_to_linear(rx, ry, 0, T_);
                        clist_[cidx].adj_cell_ids_[adj_cell_id] = cjdx;
                        adj_cell_id++;
                    }
                }
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "AdjCellList::CreateSubstructure: " << std::chrono::duration<double, std::milli> (end-start).count() << "ms\n";
    printf("AdjCellList::CreateSubstructure end\n");
}


// Update the cell list 
void
AdjCellList::UpdateCellList() {
    int midx = 0;

    // Clear the numbers in the cells
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        clist_[cidx].nparticles_ = 0;
    }

    // Now update it
    for (int i = 0; i < nparticles_; ++i) {
        int idx;
        int cx[3] = {0, 0, 0};
        const double origin[3] = {0.0, 0.0, 0.0};
        const double origins[3] = {0.0, 0.0, 0.0};
        double dr[3] = {0.0, 0.0, 0.0};
        double dr2[3] = {0.0, 0.0, 0.0};

        auto part = simples_[i];
        auto pr = part->GetPosition();
        auto prs = part->GetScaledPosition();
        // XXX: CJE must compute for arbitrary particle, maybe multiple cells (Rod)
        // Get the periodic boundary conditions on this particle
        min_distance_point_point(ndim_, nperiodic_, space_->unit_cell,
                origin, origins, pr, prs, dr, dr2);
        for (int idim = 0; idim < ndim_; ++idim) {
            cx[idim] = floor((dr[idim] + boxby2_[idim])/S_[idim]);
        }
        int cidx = cytohelpers::cell_vec_to_linear(cx[0],cx[1],cx[2], T_);
        //printf("p{%d:%d}(%2.2f, %2.2f, %2.2f) -> cell{%d}(%d,%d,%d)\n",
        //        i, oidx, part->GetPosition()[0], part->GetPosition()[1], part->GetPosition()[2],
        //        cidx, cx[0], cx[1], cx[2]);
        idx = clist_[cidx].nparticles_;
        clist_[cidx].idxlist_[idx] = i; // store the location of the particle!!!!!
        p_c_[i] = cidx; // Set the particle id->cell id
        ++idx;
        clist_[cidx].nparticles_ = idx;
        if (idx > midx) midx = idx;
    }

    if (midx > nidx_) {
        printf("Overflow in cell list: %d/%d particles/cells\n", midx, nidx_);
        exit(1);
    }
}


// print
void
AdjCellList::print() {
    printf("********\n");
    printf("%s ->\n", name_.c_str());
    printf("\t%dx%dx%d=%d cells of lengths {%.2f, %.2f, %.2f} "
           "with %d adj. cells, and %d particles/cell (n: %d).\n", T_[0], T_[1], T_[2], ncells_, 
           S_[0], S_[1], S_[2], nadj_, nidx_, nparticles_);
    printf("\t{rcut: %2.2f}, {skin: %2.2f} = {rbuff: %2.2f}\n", rcut_, skin_, rbuff_);
}


// dump gory details
void
AdjCellList::dump() {
    #ifdef DEBUG
    printf("********\n");
    printf("%s -> dump\n", name_.c_str());
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        int cx[3] = {0, 0, 0};
        cytohelpers::cell_linear_to_vec(cidx, T_, cx);
        auto cell1 = clist_[cidx];
        printf("cell{%d}[%d,%d,%d](n: %d) -> [", cidx, cx[0], cx[1], cx[2], cell1.nparticles_);
        for (int idx = 0; idx < cell1.nparticles_; ++idx) {
            printf("%d,", cell1.idxlist_[idx]);
        }
        printf("]\n");
    }
    #endif
}
