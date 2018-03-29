#include "ix_engine.h"

IxEngine::IxEngine(std::vector<Object> * o) : pt(o,&nlist) {
  objs = o;
}

void IxEngine::Init() {
  pt.AllocateCellList();
  pt.AssignCells();
  pt.CreatePairsCellList();
}

bool IxEngine::CheckUpdate() {
  double dr_max = 0;
#pragma omp parallel
  {
    double dr_local = 0;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      double dr = (*objs)[i].GetDr();
      if (dr > dr_local) dr_local = dr;
    }
#pragma omp critical
    if (dr_local > dr_max) dr_max = dr_local;
  }
  if (dr_max > 0.5*RCUT) {
    return true;
  }
  return false;
}

void IxEngine::UpdatePos0s() {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      (*objs)[i].UpdatePos0();
    }
  }
}

void IxEngine::TrackParticles() {
  if (CheckUpdate()) {
    UpdatePos0s();
    pt.AssignCells();
    pt.CreatePairsCellList();
  }
}

void IxEngine::InteractPairs() {
  int total=0;
  int nix = nlist.size();
  if (DEBUG) printf("Created %d potential interaction pairs\n",nix);
#pragma omp parallel
  {
    int count = 0;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = nix*tid/n_threads;
    int high = nix*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      Interact(nlist[i].first,nlist[i].second,&count);
    }
#pragma omp critical
    total += count;
  }
  if (DEBUG) printf("Counted %d true interactions\n",total);
}

void IxEngine::Cleanup() {
  pt.DeallocateCellList();
}

void IxEngine::Interact(int i, int j, int *count) { 
  double dx = ((*objs)[i].spos[0] - (*objs)[j].spos[0]);
  double dy = ((*objs)[i].spos[1] - (*objs)[j].spos[1]);
  double f_mag = sqrt(dx*dx+dy*dy);
  if (f_mag < 1) {
    //printf("{%d,%d}: %2.2f\n",i,j,f_mag);
    //std::lock_guard<std::mutex> lk(sim_mtx);
    (*count)++;
  }
}

