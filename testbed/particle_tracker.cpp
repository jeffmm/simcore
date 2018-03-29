#include "particle_tracker.h"

ParticleTracker::ParticleTracker(std::vector<Object> * o, std::vector<ix_pair> * n) {
  objs = o;
  nlist = n;
}

void ParticleTracker::AllocateCellList() {
  if (DEBUG) printf("Allocating cell lists\n");
  cell_length_1d = RCUT;
  n_cells_1d = (int) floor (BOX_SIZE/cell_length_1d);
  cell_length_1d = (double) BOX_SIZE/n_cells_1d;
  printf("Cell length: %2.2f\n",cell_length_1d);
  int third_dim = (NDIM==3 ? n_cells_1d : 1);
  clist = new Cell ** [n_cells_1d];
  for (int i=0;i<n_cells_1d;++i) {
    clist[i] = new Cell * [n_cells_1d];
    for (int j=0;j<n_cells_1d;++j) {
      clist[i][j] = new Cell [third_dim];
    }
  }
}

void ParticleTracker::DeallocateCellList() {
  for (int i=0;i<n_cells_1d;++i) {
    for (int j=0;j<n_cells_1d;++j) {
      delete[] clist[i][j];
    }
    delete[] clist[i];
  }
  delete[] clist;
}

void ParticleTracker::ClearCells() {
  if (DEBUG) printf("Clearing cells\n");
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_cells_1d*tid/n_threads;
    int high = n_cells_1d*(tid+1)/n_threads;
    int khigh = (NDIM==3 ? n_cells_1d : 1);
    for (int i=low; i<high; ++i) {
      for (int j=0;j<n_cells_1d;++j) {
        for (int k=0; k<khigh; ++k) {
          clist[i][j][k].objs.clear();
        }
      }
    }
  }
}

void ParticleTracker::AssignCells() {
  ClearCells();
  if (DEBUG) printf("Assigning cells\n");
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      double x = (*objs)[i].spos[0] + 0.5;
      //double x = (*objs)[i].spos[0] + 0.5*BOX_SIZE;
      int r = (int) floor (n_cells_1d*x);
      //int r = (int) floor (x/cell_length_1d);
      double y = (*objs)[i].spos[1] + 0.5;
      //double y = (*objs)[i].spos[1] + 0.5*BOX_SIZE;
      int s = (int) floor(n_cells_1d*y);
      //int s = (int) floor(y/cell_length_1d);
      int t = 0;
      if (NDIM == 3) {
        //double z = (*objs)[i].spos[2] + 0.5*BOX_SIZE;
        double z = (*objs)[i].spos[2] + 0.5;
        t = (int) floor (n_cells_1d*z);
        //t = (int) floor (z/cell_length_1d);
      }
      clist[r][s][t].AddObj(i);
    }
  }
}

void ParticleTracker::CreatePairsCellList() {
  if (DEBUG) printf("Creating pairs\n");
  nlist->clear();
#pragma omp parallel
  {
    std::vector<ix_pair> nlist_local;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_cells_1d*tid/n_threads;
    int high = n_cells_1d*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      for (int iprime = i; iprime < i+2; ++iprime) {
        int ii = iprime;
        if (ii == n_cells_1d && NPER > 0) ii = 0;
        else if (ii == n_cells_1d) continue;
        for (int j=0;j<n_cells_1d;++j) {
          for (int jprime = j-1; jprime < j+2; ++jprime) {
            int jj = jprime;
            /* We want to skip the cell below us, as that pairing will
               be taken care of by that cell, but we want to include
               the cell below and to the right */
            if (iprime != i || jprime != j-1) {
              /* Wrap for periodic boundary conditions */
              if (jj == n_cells_1d  && NPER > 1) jj = 0;
              else if (jj == n_cells_1d) continue;
              if (jj < 0 && NPER > 1) jj = n_cells_1d-1;
              else if (jj < 0) continue;
              if (NDIM<3) {
                clist[i][j][0].PairObjs(&(clist[ii][jj][0]),&nlist_local);
              }
            }
            if (NDIM<3) continue;
            for (int k=0; k<n_cells_1d; ++k) {
              for (int kprime=k; kprime<k+2; ++kprime) {
                int kk = kprime;
            /* We want to skip the cells for xprime<x, in the current y and z 
               as those pairing will be taken care of by those cells, but we
               want to include the cells xprime<x for yprime>y and zprime>z */
                if (iprime != i || kprime != k-1 || jprime != j-1) {
                  /* Wrap for periodic boundary conditions */
                  if (kk == n_cells_1d && NPER > 2) kk = 0;
                  else if (kk == n_cells_1d) continue;
                  if (kk < 0 && NPER > 2) kk = n_cells_1d-1;
                  else if (kk < 0) continue;
                  clist[i][j][k].PairObjs(&(clist[ii][jj][kk]),&nlist_local);
                }
              }
            }
          }
        }
      }
    }
#pragma omp critical
    nlist->insert(nlist->end(), nlist_local.begin(), nlist_local.end());
  }
}

