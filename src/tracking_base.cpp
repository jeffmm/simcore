// Implementation of tracking base class

#include "tracking_base.h"

// Initialize
void
TrackingBase::Init(space_struct* pSpace,
                   std::vector<SpeciesBase*> *pSpecies,
                   std::vector<Simple*> *pSimples,
                   PotentialManager *pPotentials,
                   double pSkin) {
    space_ = pSpace;
    simples_ = pSimples;
    species_ = pSpecies;
    ndim_ = space_->n_dim;
    nperiodic_ = space_->n_periodic;
    skin_ = pSkin;
    potentials_ = pPotentials;
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
    nsimples_ = (int)simples_->size();
    nupdates_ = 0;

    rid_interactions_ = new std::unordered_set<std::pair<int, int>, hashh::pair_hash>[nthreads_+1];
    rid_self_check_ = new std::vector<bool>();
    //unique_rids_ = new std::unordered_set<int>();
    unique_rids_ = new std::set<int>();
    rid_check_local_ = new std::unordered_set<int>*[nthreads_];
    for (int ithread = 0; ithread < nthreads_; ++ithread) {
      rid_check_local_[ithread] = new std::unordered_set<int>();
    }
}

void TrackingBase::Rebuild(nl_list **pNeighbors) {
  neighbors_ = (*pNeighbors);
  nsimples_ = (int)simples_->size();
}
