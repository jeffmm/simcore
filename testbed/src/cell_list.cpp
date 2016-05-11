#include "cell_list.h"

CellList::CellList() {}

// First, generate the cells
// Independent of particle type, just need the best
// rcutoff to generate this
void
CellList::generateCells(int npart, double box, double rcut) {
    int i, j, k;
    int npair_loc;
    double boxby2, boxoffs;
    boxby2 = 0.5 * box;
    nparticles = npart;

    ngrid = floor(cellrat * box / rcut);
    ncell = ngrid*ngrid*ngrid;
    delta = box / ngrid;
    boxoffs = boxby2 - 0.5*delta;

    // Set storage
    clist.clear();
    plist.clear();
    clist.resize(ncell);
    plist.resize(2*ncell*ncell);

    // Allocate index lists within cells, cell density < 2x avg. density
    nidx = 2*nparticles/ncell + 2;
    nidx = ((nidx/2)+1) * 2;
    for (int i = 0; i < ncell; ++i) {
        clist[i].idxlist.clear();
        clist[i].ptidxlist.clear();
        clist[i].idxlist.resize(nidx);
        clist[i].ptidxlist.resize(nidx);
    }

    // Build the cell pair list
    npair_loc = 0;
    for (i = 0; i < ncell-1; ++i) {
        double x1, x2, y1, y2, z1, z2, rx, ry, rz;

        k = i/ngrid/ngrid;
        x1 = k*delta - boxoffs;
        y1 = ((i-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
        z1 = (i % ngrid)*delta - boxoffs;

        for (j = i+1; j < ncell; ++j) {
            k = j/ngrid/ngrid;
            x2 = k*delta - boxoffs;
            y2 = ((j-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
            z2 = (j % ngrid)*delta - boxoffs;

            rx = buffmd::pbc(x1 - x2, boxby2, box);
            ry = buffmd::pbc(y1 - y2, boxby2, box);
            rz = buffmd::pbc(z1 - z2, boxby2, box);

            // Check the cells on a line that are too far apart
            if (fabs(rx) > rcut + delta) continue;
            if (fabs(ry) > rcut + delta) continue;
            if (fabs(rz) > rcut + delta) continue;

            // Check for cells in a plane that are too far apart
            if (sqrt(rx*rx + ry*ry) > (rcut + sqrt(2.0)*delta)) continue;
            if (sqrt(rx*rx + rz*rz) > (rcut + sqrt(2.0)*delta)) continue;
            if (sqrt(ry*ry + rz*rz) > (rcut + sqrt(2.0)*delta)) continue;

            // Other cells too far apart in 3d
            if (sqrt(rx*rx + ry*ry + rz*rz) > (sqrt(3.0)*delta + rcut)) continue;

            // Cells are close enough, add
            plist[2*npair_loc    ] = i;
            plist[2*npair_loc + 1] = j;
            ++npair_loc;
        } // j
    } // i
    npair = npair_loc;
    printf("Cell list has %dx%dx%d=%d cells with %d/%d pairs and "
            "%d particles/celllist.\n", ngrid, ngrid, ngrid, ncell,
            npair, ncell*(ncell-1)/2, nidx);
}

// Update the cell list based on base_particles 
void
CellList::updateCells(double box, std::vector<particle*>* particles) {
    // Reset the cell list
    int k, m, n, j;

    double boxby2 = 0.5*box;
    for (int i = 0; i < ncell; ++i) {
        clist[i].nparticles = 0;
    }

    midx = 0;
    for(int i = 0; i < nparticles; ++i) {
        int idx;

        auto p = (*particles)[i];
        k = floor((buffmd::pbc(p->x[0], boxby2, box) + boxby2)/delta);
        m = floor((buffmd::pbc(p->x[1], boxby2, box) + boxby2)/delta);
        n = floor((buffmd::pbc(p->x[2], boxby2, box) + boxby2)/delta);
        j = ngrid*ngrid*k + ngrid*m + n;

        idx = clist[j].nparticles;
        clist[j].idxlist[idx] = i;
        clist[j].ptidxlist[idx] = p->sid;
        ++idx;
        clist[j].nparticles = idx;
        if (idx > midx) midx = idx;
    }

    if (midx > nidx) {
        printf("Overflow in cell list: %d/%d particles/cells\n", midx, nidx);
        exit(1);
    }
}


// Return the number of cells
int
CellList::nCells() {
    return ncell;
}


// Return the number of pairs
int
CellList::nPairs() {
    return npair;
}


// Get a particular cell
cell_t*
CellList::getCell(int idx) {
    return &clist[idx];
}


// Get the pair cell
int
CellList::getPair(int idx) {
    return plist[idx];
}
