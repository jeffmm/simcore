// Implementation of tracking base class

#include "tracking_base.h"

// Initialize
void
TrackingBase::Init(space_struct* pSpace, std::vector<Simple*> *pSimples, double pSkin) {
    space_ = pSpace;
    simples_ = pSimples;
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
    nsimples_ = (int)simples_->size();
}

