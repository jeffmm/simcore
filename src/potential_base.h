#ifndef _CYTOSCORE_POTENTIAL_BASE_H_
#define _CYTOSCORE_POTENTIAL_BASE_H_

#include "auxiliary.h"

// Potential base class for all external(or maybe even internal) potentials
class PotentialBase {
  protected:
    int ndim_;
    double rcut_, rcut2_; // Cutoff radius
    space_struct* space_;
  public:
    PotentialBase(space_struct* pSpace, double pRcut) : space_(pSpace), rcut_(pRcut) {
        rcut2_ = rcut_ * rcut_;
        ndim_ = space_->n_dim;
    }
    virtual ~PotentialBase() {}
    virtual void CalcPotential(double* x,
                               double* xs,
                               double* y,
                               double* ys,
                               double* fpote) {}
    virtual void Print() {
        printf("\t{rcut:%2.2f}\n", rcut_);
        for (int i = 0; i < ndim_; ++i) {
            printf("\t");
            for (int j = 0; j < ndim_; ++j) {
                printf("%2.2f\t", space_->unit_cell[i][j]);
            }
            printf("\n");
        }
    }
};

#endif
