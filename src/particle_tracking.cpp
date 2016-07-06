// Implementation file for particle tracking (bookkeeping)

#include <climits>

#include "particle_tracking.h"

// Pass in the main system parameters
void ParticleTracking::Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin, FTYPE pFtype) {
  space_ = pSpace;
  species_ = pSpecies;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  skin_ = pSkin;
  ftype_ = pFtype;
  for (int i = 0; i < ndim_; ++i) {
    box_[i] = space_->unit_cell[i][i];
  }
  printf("********\n");
  printf("Initializing Particle Tracking Scheme ->\n");

  for (int i = 0; i < ndim_; ++i) {
    box_[i] = space_->unit_cell[i][i];
  }

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  {
    if (0 == omp_get_thread_num()) {
      nthreads_ = omp_get_num_threads();
    }
  }
  #else
  nthreads_ = 1;
  #endif

  // Figure out what kind of substructure we're using
  switch (ftype_) {
    case FTYPE::allpairs:
      printf("\tUsing all pairs particle tracking substructure\n");
      tracking_ = trackingFactory<TrackingAllPairs>();
      break;
    case FTYPE::neighborallpairs:
      printf("\tUsing neighbor list-all pairs tracking substructure\n");
      tracking_ = trackingFactory<TrackingNeighborListAP>();
      break;
    default:
      printf("Must specify a particle tracking substructure, exiting!\n");
      exit(1);
  }

  // Load the simples directly
  LoadSimples();
}

// Load our simples array
void ParticleTracking::LoadSimples() {
  nsys_ = species_->size();
  simples_.clear();
  for (auto it = species_->begin(); it != species_->end(); ++it) {
    std::vector<Simple*> sim_vec = (*it)->GetSimples();
    simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
  }
  nsimples_ = (int)simples_.size();
}

// Attach to the potential manager
void ParticleTracking::InitPotentials(PotentialManager* pPotentials) {
  potentials_ = pPotentials;
  rcut_ = potentials_->GetMaxRCut();
}

// Initialize the tracking
void ParticleTracking::InitTracking() {
  neighbors_ = new nl_list[nsimples_];

  // Create the neighbor list based on the tracking criteria
  // IE Brute
  tracking_->Init(space_, species_, &simples_, skin_);
  tracking_->CreateSubstructure(rcut_, &neighbors_);
  tracking_->UpdateTracking(true);
}

// Check if we need an update
// and do it if needed
void ParticleTracking::UpdateTracking(bool pForceUpdate) {
  tracking_->UpdateTracking(pForceUpdate);
}

// Check the overlaps
void ParticleTracking::CheckOverlaps(int pMaxOverlaps) {
  // XXX CJE fix this to use substructure
  // very yucky brute force method of checking overlaps
  // load all the particles
  bool overlap = true;
  int num = 0;
  do {
    num++;
    overlap = false;
    LoadSimples();
    for (int idx = 0; idx < simples_.size() - 1 && !overlap; ++idx) {
      for (int jdx = idx + 1; jdx < simples_.size() && !overlap; ++jdx) {
        auto part1 = simples_[idx];
        auto part2 = simples_[jdx];

        if (part1->GetCID() == part2->GetCID()) continue;

        // Check to see if interacting
        PotentialBase *pot = potentials_->GetPotential(part1->GetSID(), part2->GetSID());
        if (pot == nullptr) continue;

        // Check to see if overlap on this is allowed (KMC stuff)
        if (pot->CanOverlap()) continue;

        // Can possibly overlap, check
        interactionmindist idm;
        MinimumDistance(part1, part2, idm, ndim_, nperiodic_, space_);
        if (idm.dr_mag2 < idm.buffer_mag2) {
          if (debug_trace) {
            printf("Overlap detected [oid: %d,%d], [idx: %d, %d], [sid: %d, %d], [cid: %d, %d] -> (%2.2f, %2.2f)\n",
                    part1->GetOID(), part2->GetOID(), idx, jdx, part1->GetSID(), part2->GetSID(), 
                    part1->GetCID(), part2->GetCID(), idm.buffer_mag2, idm.dr_mag2);
          }
          overlap = true;
          auto sid = part2->GetSID();
          auto cid = part2->GetCID();
          for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
            if ((*spec)->GetSID() == sid) {
              (*spec)->ReInit(cid);
            }
          }
        }
      }
    }
    if (num > pMaxOverlaps)
        error_exit("ERROR: Too many overlaps detected. Check packing ratio for objects.\n");
  } while (overlap);
}

// Print out interesting information
void ParticleTracking::Print() {
  printf("********\n");
  printf("%s ->\n", tracking_->Name().c_str());
  printf("\t{rcut: %2.2f}, {skin: %2.2f}\n", rcut_, skin_);
  int ntotlist = 0, maxlist = 0, minlist = INT_MAX;
  for (int i = 0; i < nsimples_ -1; ++i) {
    ntotlist += neighbors_[i].size();
    maxlist = std::max(maxlist, (int)neighbors_[i].size());
    minlist = std::min(minlist, (int)neighbors_[i].size());
  }
  printf("\tStats: {min: %d}, {max: %d}, {avg: %2.2f}\n", minlist, maxlist, (float)ntotlist/(float)nsimples_);
  tracking_->print();
}

// Dump the gory details of the neighbor list
void ParticleTracking::Dump() {
  #ifdef DEBUG
  printf("********\n");
  printf("%s -> dump\n", tracking_->Name().c_str());
  for (int idx = 0; idx < nsimples_; ++idx) {
    printf("\t[%d] -> [", idx);
    for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
      int jdx = nldx->idx_;
      printf("%d,", jdx);
    }
    printf("]\n");
  }
  tracking_->dump();
  #endif
}
