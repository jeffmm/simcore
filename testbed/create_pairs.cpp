#include "simulation.h"

// First attempt, pair via box
void Simulation::CreatePairs() {
  //int n_cell_1d = (int) floor(BOX_SIZE/LCELL);
  //double cell_length = BOX_SIZE/n_cell_1d;
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      double x = objs[i].pos[0];
      double y = objs[i].pos[1];
      for (int j=i+1; j<NOBJS; ++j) {
        if (abs(objs[j].pos[0] - x) < RCUT && abs(objs[j].pos[1] - y) < RCUT) {
          objs[i].nlist.push_back(j);
        }
      }
    }
  }
}

void Simulation::PrintPairs() {
  for (int i=0; i<NOBJS; ++i) {
    printf("  %d interacts with { ",i);
    for (auto it=objs[i].nlist.begin(); it!=objs[i].nlist.end(); ++it) {
      printf("%d ",(*it));
    }
    printf("}\n");
  }
}
