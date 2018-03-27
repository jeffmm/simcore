#include "fff_simulation.h"

void FFFSim::Run() {
  std::cout << "Running simulation on " << omp_get_max_threads() << " threads" <<  std::endl;
  if (DEBUG) printf("Creating Objects\n");
  objs.resize(NOBJS);
  if (DEBUG) printf("Initiating RNG\n");
  InitRNG(); // Init RNGs
  if (DEBUG) printf("Initializing positions\n");
  InitPositions(); // Init object positions
  ix_engine.Init();
  for (int i=0;i<1000;++i) {
    UpdatePositions();
    ix_engine.TrackParticles();
    ix_engine.InteractPairs();
  }
  ix_engine.Cleanup();
  printf("Done\n");
  //PrintPairs();
}

void FFFSim::UpdatePositions() {
  if (DEBUG) printf("Updating positions\n");
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      objs[i].UpdatePosition();
    }
  }
}

// Initialize random number generators
void FFFSim::InitRNG() {
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

void FFFSim::InitPositions() {
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    int low = NOBJS*tid/n_threads;
    int high = NOBJS*(tid+1)/n_threads;
    for (int i=low; i<high; ++i) {
      objs[i].InitRandomPosition();;
    }
  }
}

// Second attempt, pair via cell lists

// First attempt, pair via box
//void FFFSim::CreatePairs() {
  //std::cout << " Running on " << omp_get_max_threads() << " threads\n";
//#pragma omp parallel
  //{
    //std::vector<ix_pair> nlist_local;
    //int tid = omp_get_thread_num();
    //int n_threads = omp_get_num_threads();
    //int low = NOBJS*tid/n_threads;
    //int high = NOBJS*(tid+1)/n_threads;
////#pragma omp for nowait //fill local vectors in parallel
    //for (int i=low; i<high; ++i) {
      //double x = objs[i].pos[0];
      //double y = objs[i].pos[1];
      //for (int j=i+1; j<NOBJS; ++j) {
        //if (abs(objs[j].pos[0] - x) < 3*RCUT && abs(objs[j].pos[1] - y) < 3*RCUT) {
          //nlist_local.push_back(std::make_pair(i,j));
        //}
      //}
    //}
//#pragma omp critical
  ////std::lock_guard<std::mutex> lk(sim_mtx);
  //nlist.insert(nlist.end(), nlist_local.begin(), nlist_local.end());
  //}
//}


