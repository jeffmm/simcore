#include "particle_tracker.h"

void ParticleTracker::Init(system_parameters * p, std::vector<Object*> * o, std::vector<ix_pair> * n) {
  params_ = p;
  objs_ = o;
  nlist_ = n;
  n_dim_ = params_->n_dim;
  n_per_ = params_->n_periodic;
  cell_length_1d_ = params_->cell_length;
  n_cells_1d_ = (int) floor (2*params_->system_radius/cell_length_1d_);
  cell_length_1d_ = (double) 2*params_->system_radius/n_cells_1d_;
  AllocateCellList();
}

void ParticleTracker::AllocateCellList() {
  DPRINTF("Allocating cell lists\n");
  printf("Cell length: %2.2f\n",cell_length_1d_);
  int third_dim = (n_dim_==3 ? n_cells_1d_ : 1);
  clist_ = new Cell ** [n_cells_1d_];
  for (int i=0;i<n_cells_1d_;++i) {
    clist_[i] = new Cell * [n_cells_1d_];
    for (int j=0;j<n_cells_1d_;++j) {
      clist_[i][j] = new Cell [third_dim];
    }
  }
}

void ParticleTracker::DeallocateCellList() {
  for (int i=0;i<n_cells_1d_;++i) {
    for (int j=0;j<n_cells_1d_;++j) {
      delete[] clist_[i][j];
    }
    delete[] clist_[i];
  }
  delete[] clist_;
}

void ParticleTracker::ClearCells() {
  DPRINTF("Clearing cells\n");
#ifdef ENABLE_OPENMP
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_cells_1d_*tid/n_threads;
    int high = n_cells_1d_*(tid+1)/n_threads;
    int khigh = (n_dim_==3 ? n_cells_1d_ : 1);
    for (int i=low; i<high; ++i) {
      for (int j=0;j<n_cells_1d_;++j) {
        for (int k=0; k<khigh; ++k) {
          clist_[i][j][k].objs_.clear();
        }
      }
    }
  }
#else
  int khigh = (n_dim_==3 ? n_cells_1d_ : 1);
  for (int i=0; i<n_cells_1d_; ++i) {
    for (int j=0;j<n_cells_1d_;++j) {
      for (int k=0; k<khigh; ++k) {
        clist_[i][j][k].objs_.clear();
      }
    }
  }
#endif
}

void ParticleTracker::AssignCells() {
  ClearCells();
  int n_objs_ = objs_->size();
  DPRINTF("Assigning cells\n");
#ifdef ENABLE_OPENMP
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_objs_*tid/n_threads;
    int high = n_objs_*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      double const * const spos = (*objs_)[i]->GetScaledPosition();
      double x = spos[0] + 0.5;
      int r = (int) floor (n_cells_1d_*x);
      double y = spos[1] + 0.5;
      int s = (int) floor(n_cells_1d_*y);
      int t = 0;
      if (n_dim_ == 3) {
        double z = spos[2] + 0.5;
        t = (int) floor (n_cells_1d_*z);
      }
      clist_[r][s][t].AddObj(i);
    }
  }
#else
  for (int i=0; i<n_objs_; ++i) {
    double const * const spos = (*objs_)[i]->GetScaledPosition();
    double x = spos[0] + 0.5;
    int r = (int) floor (n_cells_1d_*x);
    double y = spos[1] + 0.5;
    int s = (int) floor(n_cells_1d_*y);
    int t = 0;
    if (n_dim_ == 3) {
      double z = spos[2] + 0.5;
      t = (int) floor (n_cells_1d_*z);
    }
    clist_[r][s][t].AddObj(i);
  }
#endif
}

void ParticleTracker::CreatePairsCellList() {
  DPRINTF("Creating pairs\n");
  nlist_->clear();
#ifdef ENABLE_OPENMP
#pragma omp parallel
  {
    std::vector<ix_pair> nlist_local;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_cells_1d_*tid/n_threads;
    int high = n_cells_1d_*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      for (int iprime = i; iprime < i+2; ++iprime) {
        int ii = iprime;
        if (ii == n_cells_1d_ && n_per_ > 0) ii = 0;
        else if (ii == n_cells_1d_) continue;
        for (int j=0;j<n_cells_1d_;++j) {
          for (int jprime = j-1; jprime < j+2; ++jprime) {
            int jj = jprime;
            /* We want to skip the cell below us, as that pairing will
               be taken care of by that cell, but we want to include
               the cell below and to the right */
            if (iprime != i || jprime != j-1) {
              /* Wrap for periodic boundary conditions */
              if (jj == n_cells_1d_  && n_per_ > 1) jj = 0;
              else if (jj == n_cells_1d_) continue;
              if (jj < 0 && n_per_ > 1) jj = n_cells_1d_-1;
              else if (jj < 0) continue;
              if (n_dim_<3) {
                clist_[i][j][0].PairObjs(&(clist_[ii][jj][0]),&nlist_local);
              }
            }
            if (n_dim_<3) continue;
            for (int k=0; k<n_cells_1d_; ++k) {
              for (int kprime=k; kprime<k+2; ++kprime) {
                int kk = kprime;
            /* We want to skip the cells for xprime<x, in the current y and z 
               as those pairing will be taken care of by those cells, but we
               want to include the cells xprime<x for yprime>y and zprime>z */
                if (iprime != i || kprime != k-1 || jprime != j-1) {
                  /* Wrap for periodic boundary conditions */
                  if (kk == n_cells_1d_ && n_per_ > 2) kk = 0;
                  else if (kk == n_cells_1d_) continue;
                  if (kk < 0 && n_per_ > 2) kk = n_cells_1d_-1;
                  else if (kk < 0) continue;
                  clist_[i][j][k].PairObjs(&(clist_[ii][jj][kk]),&nlist_local);
                }
              }
            }
          }
        }
      }
    }
#pragma omp critical
    nlist_->insert(nlist_->end(), nlist_local.begin(), nlist_local.end());
  }
#else
  for (int i=0; i<n_cells_1d_; ++i) {
    for (int iprime = i; iprime < i+2; ++iprime) {
      int ii = iprime;
      if (ii == n_cells_1d_ && n_per_ > 0) ii = 0;
      else if (ii == n_cells_1d_) continue;
      for (int j=0;j<n_cells_1d_;++j) {
        for (int jprime = j-1; jprime < j+2; ++jprime) {
          int jj = jprime;
          /* We want to skip the cell below us, as that pairing will
             be taken care of by that cell, but we want to include
             the cell below and to the right */
          if (iprime != i || jprime != j-1) {
            /* Wrap for periodic boundary conditions */
            if (jj == n_cells_1d_  && n_per_ > 1) jj = 0;
            else if (jj == n_cells_1d_) continue;
            if (jj < 0 && n_per_ > 1) jj = n_cells_1d_-1;
            else if (jj < 0) continue;
            if (n_dim_<3) {
              clist_[i][j][0].PairObjs(&(clist_[ii][jj][0]),nlist_);
            }
          }
          if (n_dim_<3) continue;
          for (int k=0; k<n_cells_1d_; ++k) {
            for (int kprime=k; kprime<k+2; ++kprime) {
              int kk = kprime;
          /* We want to skip the cells for xprime<x, in the current y and z 
             as those pairing will be taken care of by those cells, but we
             want to include the cells xprime<x for yprime>y and zprime>z */
              if (iprime != i || kprime != k-1 || jprime != j-1) {
                /* Wrap for periodic boundary conditions */
                if (kk == n_cells_1d_ && n_per_ > 2) kk = 0;
                else if (kk == n_cells_1d_) continue;
                if (kk < 0 && n_per_ > 2) kk = n_cells_1d_-1;
                else if (kk < 0) continue;
                clist_[i][j][k].PairObjs(&(clist_[ii][jj][kk]),nlist_);
              }
            }
          }
        }
      }
    }
  }
#endif
}

void ParticleTracker::Check

void ParticleTracker::Clear() {
  DeallocateCellList();
}

double ParticleTracker::GetCellLength() {
  return cell_length_1d_;
}

void ParticleTracker::CreatePartialPairsCellList() {

}
