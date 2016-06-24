// Implementation of forcesubstructure base class

#include "force_substructure_base.h"

// Initialize
void
ForceSubstructureBase::Init(space_struct* pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin) {
    space_ = pSpace;
    species_ = pSpecies;
    ndim_ = space_->n_dim;
    nperiodic_ = space_->n_periodic;
    skin_ = pSkin;
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
}

// Load the simple particles into the master vector
void
ForceSubstructureBase::LoadFlatSimples(std::vector<Simple*> pSimples) {
 
    // Load everything from the our given simples
    simples_.clear();
    simples_.insert(simples_.end(), pSimples.begin(), pSimples.end());
    nparticles_ = (int)simples_.size();
}
