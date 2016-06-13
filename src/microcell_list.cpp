// implementation for microcell list

#include "microcell_list.h"

#include "minimum_distance.h"

void
MicrocellList::CreateSubstructure(double pRcut) {
    double boxoffs[3];

    printf("MicrocellList::CreateSubstructure\n");
    rcut_ = pRcut;
    rbuff_ = rcut_ + skin_;
    p_c_.clear();
    p_c_.resize(nparticles_);

    // Compute the number of cells on a side
    for (int i = 0; i < ndim_; ++i) {
        boxby2_[i] = 0.5*box_[i];
        T_[i] = floor(cellrat_ * box_[i] / rbuff_);
        if (T_[i] < 3) {
            T_[i] = 3;
        }
    }

    // Compute the side length
    for (int i = 0; i < ndim_; ++i) {
        S_[i] = box_[i] / T_[i];
        boxoffs[i] = boxby2_[i] - 0.5*S_[i];
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
    plist_.clear();
    plist_.resize(2*ncells_*ncells_);

    // Allocate the list index within cells
    for (int i = 0; i < ncells_; ++i) {
        clist_[i].cell_id_ = i;
        clist_[i].idxlist_.clear();
        clist_[i].idxlist_.resize(nidx_);
    }

    // Build the cell pair list finally
    npairs_ = 0;
    for (int cidx = 0; cidx < ncells_ -1; ++cidx) {
        int cx1[3] = {0, 0, 0};
        double x1[3] = {0.0, 0.0, 0.0};
        double x2[3] = {0.0, 0.0, 0.0};

        // First cell
        cytohelpers::cell_linear_to_vec(cidx, T_, cx1);
        for (int idim = 0; idim < ndim_; ++idim) {
            x1[idim] = cx1[idim]*S_[idim] - boxoffs[idim];
        }

        // Second cell
        for (int cjdx = cidx+1; cjdx < ncells_; ++cjdx) {
            int cx2[3] = {0, 0, 0};
            double dr[3] = {0.0, 0.0, 0.0};
            double ds[3] = {0.0, 0.0, 0.0};

            cytohelpers::cell_linear_to_vec(cjdx, T_, cx2);
            for (int idim = 0; idim < ndim_; ++idim) {
                x2[idim] = cx2[idim]*S_[idim] - boxoffs[idim];
                dr[idim] = x1[idim] - x2[idim];
            }

            // Compute the separation *with* periodic boundary conditions
            periodic_boundary_conditions(nperiodic_, space_->unit_cell, space_->unit_cell_inv,
                                         dr, ds);

            // Check if the cells are too far apart in all ways
            // 1d
            if (fabs(dr[0]) > rbuff_ + S_[0]) continue;
            if (fabs(dr[1]) > rbuff_ + S_[1]) continue;
            if (fabs(dr[2]) > rbuff_ + S_[2]) continue;

            // 2d
            if (sqrt(dr[0]*dr[0] + dr[1]*dr[1]) > (rbuff_ + sqrt(S_[0]*S_[0] + S_[1]*S_[1]))) continue;
            if (sqrt(dr[0]*dr[0] + dr[2]*dr[2]) > (rbuff_ + sqrt(S_[0]*S_[0] + S_[2]*S_[2]))) continue;
            if (sqrt(dr[1]*dr[1] + dr[2]*dr[2]) > (rbuff_ + sqrt(S_[1]*S_[1] + S_[2]*S_[2]))) continue;

            // 3d
            if (sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) > (sqrt(S_[0]*S_[0] + S_[1]*S_[1] + S_[2]*S_[2]) + rbuff_)) continue;
            // Add to the pair list if close enough
            plist_[2*npairs_  ] = cidx;
            plist_[2*npairs_+1] = cjdx;
            ++npairs_;
        }
    }

    #ifdef DEBUG
    // Write out the cell locations and extent
    for (int cidx = 0; cidx < ncells_; ++cidx) {
        int cx[3] = {0, 0, 0};
        double x[3] = {0.0, 0.0, 0.0};
        cytohelpers::cell_linear_to_vec(cidx, T_, cx);
        for (int idim = 0; idim < ndim_; ++idim) {
            x[idim] = cx[idim]*S_[idim] - boxoffs[idim];
        }
        printf("cell{%d}[%d,%d,%d] -> (%2.2f, %2.2f, %2.2f), ",
                cidx, cx[0], cx[1], cx[2], x[0], x[1], x[2]);
        printf("lo(%2.2f, %2.2f, %2.2f) - hi(%2.2f, %2.2f, %2.2f)\n",
                x[0]-0.5*S_[0], x[1]-0.5*S_[1], x[2]-0.5*S_[2],
                x[0]+0.5*S_[0], x[1]+0.5*S_[1], x[2]+0.5*S_[2]);
    }
    #endif

    printf("********\n");
    printf("Microcell list (ndim=%d) has %dx%dx%d=%d cells of lengths {%.2f, %.2f, %.2f} "
           "with %d/%d pairs and %d particles/cell.\n", ndim_, T_[0], T_[1], T_[2], ncells_, 
           S_[0], S_[1], S_[2], npairs_, ncells_*(ncells_-1)/2, nidx_);
}


// Update the cell list 
void
MicrocellList::UpdateCellList() {
    int midx = 0;
    int cx[3];

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
        auto oidx = part->GetOID();
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
