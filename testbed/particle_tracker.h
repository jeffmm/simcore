#include "object.h"
#include "omp.h"
#include <iostream>
#include <math.h>
#include "cell_list.h"

class ParticleTracker {
  private:
    double cell_length_1d;
    int n_cells_1d;
    Cell *** clist;
    std::vector<ix_pair> * nlist;
    std::vector<Object> * objs;
  public:
    ParticleTracker(std::vector<Object> * o, std::vector<ix_pair> * n);
    void AllocateCellList();
    void DeallocateCellList();
    void CreatePairsCellList();
    void AssignCells();
    void ClearCells();
};
