// Most basic Cell List implementation I can think of
// Operates on base particles, with some given cutoff
#ifndef BUFFMD_CELLLIST_H_
#define BUFFMD_CELLLIST_H_

#include <vector>
#include <cmath>

#include "particle.h"
#include "species.h"

// Each actual cell has information about the ids of the particles
struct _cell {
    int nparticles;
    int owner;
    // Particle id list
    std::vector<int> idxlist;
    // Particle type list
    std::vector<int> ptidxlist;
};
typedef struct _cell cell_t;

class CellList {
public:
    CellList();
    void generateCells(int npart, double box, double rcut);
    void updateCells(double box, std::vector<particle*>* particles);
    int nCells();
    int nPairs();
    cell_t* getCell(int idx);
    int getPair(int idx);

protected:
    int ngrid;
    int ncell;
    int npair;
    int nidx, midx;
    int nparticles;

    double delta;
    const double cellrat = 2.0;

    std::vector<cell_t> clist;
    std::vector<int> plist;
};

#endif
