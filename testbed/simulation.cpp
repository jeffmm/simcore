#include "simulation.h"

void Simulation::Run() {
  printf("Creating Objects\n");
  objs.resize(NOBJS);
  AllocateCellList();
  printf("Initiating RNG\n");
  InitRNG(); // Init RNGs
  printf("Initializing positions\n");
  InitPositions(); // Init object positions
  for (int i=0;i<10;++i) {
    printf("Assigning cells\n");
    AssignCells();
    printf("Creating pairs\n");
    //CreatePairs();
    CreatePairsCellList();
    printf("Interacting\n");
    InteractPairs();
  }
  DeallocateCellList();
  printf("Done\n");
  //PrintPairs();
}

// Initialize random number generators
void Simulation::InitRNG() {
  gsl_rng_env_setup();
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      objs[i].rng = gsl_rng_alloc(gsl_rng_taus);
      gsl_rng_set(objs[i].rng, (i+1)*SEED);
    }
  }
}

void Simulation::AllocateCellList() {
  cell_length_1d = RCUT;
  n_cells_1d = (int) floor (BOX_SIZE/cell_length_1d);
  cell_length_1d = (double) BOX_SIZE/n_cells_1d;
  printf("Cell length: %2.2f\n",cell_length_1d);
  clist = new Cell ** [n_cells_1d];
  for (int i=0;i<n_cells_1d;++i) {
    clist[i] = new Cell * [n_cells_1d];
    for (int j=0;j<n_cells_1d;++j) {
      clist[i][j] = new Cell [NDIM==3 ? (n_cells_1d) : 1];
    }
  }
}

void Simulation::DeallocateCellList() {
  for (int i=0;i<n_cells_1d;++i) {
    for (int j=0;j<n_cells_1d;++j) {
      delete[] clist[i][j];
    }
    delete[] clist[i];
  }
  delete[] clist;
}

void Simulation::InitPositions() {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      objs[i].SetRandomPosition();;
    }
  }
}

void Simulation::AssignCells() {
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
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      double x = objs[i].pos[0] + 0.5*BOX_SIZE;
      int r = (int) floor (x/cell_length_1d);
      double y = objs[i].pos[1] + 0.5*BOX_SIZE;
      int s = (int) floor(y/cell_length_1d);
      int t = 0;
      if (NDIM == 3) {
        double z = objs[i].pos[2] + 0.5*BOX_SIZE;
        t = (int) floor (z/cell_length_1d);
      }
      clist[r][s][t].AddObj(i);
    }
  }
}

// Second attempt, pair via cell lists
void Simulation::CreatePairsCellList() {
  std::cout << " Running on " << omp_get_max_threads() << " threads\n";
  nlist.clear();
#pragma omp parallel
  {
    std::vector<ix_pair> nlist_local;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = n_cells_1d*tid/n_threads;
    int high = n_cells_1d*(tid+1)/n_threads;
//#pragma omp for nowait //fill local vectors in parallel
    for (int i=low; i<high; ++i) {
      for (int iprime = i; iprime < i+2; ++iprime) {
        if (iprime == n_cells_1d) continue;//iprime = 0;
        for (int j=0;j<n_cells_1d;++j) {
          for (int jprime = j-1; jprime < j+2; ++jprime) {
            if (iprime == i && jprime == j-1) continue;
            if (jprime == n_cells_1d ) continue;//jprime = 0;
            if (jprime < 0) continue;//jprime = n_cells_1d-1;
            //if (NDIM<3) {
              clist[i][j][0].PairObjs(&(clist[iprime][jprime][0]),&nlist_local);
              //continue;
            //}
            //for (int k=0; k<n_cells_1d; ++k) {
              //for (int kk=k; kk<k+2; ++kk) {
                //if (kk == n_cells_1d) kk = 0;
                //// First take care of pairing objects in this cell
                //clist[i][j][k].PairObjs(&(clist[ii][jj][kk]),&nlist_local);
              //}
            //}
          }
        }
      }
    }
#pragma omp critical
    //std::lock_guard<std::mutex> lk(sim_mtx);
    nlist.insert(nlist.end(), nlist_local.begin(), nlist_local.end());
  }
}

// First attempt, pair via box
void Simulation::CreatePairs() {
  std::cout << " Running on " << omp_get_max_threads() << " threads\n";
#pragma omp parallel
  {
    std::vector<ix_pair> nlist_local;
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
//#pragma omp for nowait //fill local vectors in parallel
    for (int i=low; i<high; ++i) {
      double x = objs[i].pos[0];
      double y = objs[i].pos[1];
      for (int j=i+1; j<NOBJS; ++j) {
        if (abs(objs[j].pos[0] - x) < 3*RCUT && abs(objs[j].pos[1] - y) < 3*RCUT) {
          nlist_local.push_back(std::make_pair(i,j));
        }
      }
    }
#pragma omp critical
  //std::lock_guard<std::mutex> lk(sim_mtx);
  nlist.insert(nlist.end(), nlist_local.begin(), nlist_local.end());
  }
}

void Simulation::PrintPairs() {
  //for (int i=0; i<NOBJS; ++i) {
    //printf("  %d interacts with { ",i);
    //for (auto it=objs[i].nlist.begin(); it!=objs[i].nlist.end(); ++it) {
      //printf("%d ",(*it));
    //}
    //printf("}\n");
  //}
}

void Simulation::InteractPairs() {
  int nix = nlist.size();
  printf("%d interactions\n",nix);
  //for (int i=0; i<nix; ++i) {
    //Interact(nlist[i].first,nlist[i].second);
  //}
  int total=0;
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
  printf("\nCount: %d\n",total);
}

void Simulation::Interact(int i, int j, int *count) { 
  double dx = (objs[i].pos[0] - objs[j].pos[0]);
  double dy = (objs[i].pos[1] - objs[j].pos[1]);
  double f_mag = sqrt(dx*dx+dy*dy);
  if (f_mag < 1) {
    //printf("{%d,%d}: %2.2f\n",i,j,f_mag);
    //std::lock_guard<std::mutex> lk(sim_mtx);
    (*count)++;
  }
}

