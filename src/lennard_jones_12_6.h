#ifndef _CYTOSCORE_LJ126_H_
#define _CYTOSCORE_LJ126_H_

#include <cmath>
#include <cstdio>

#include "potential_base.h"

class LJ126 : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
  public:
    LJ126(double pEps, double pSigma, space_struct* pSpace, double pRcut) : PotentialBase(pSpace, pRcut), eps_(pEps), sigma_(pSigma) {
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
    virtual void Print() {
        std::cout << "Lennard-Jones 12-6 potential:\n";
        PotentialBase::Print();
        std::cout << "\t{eps:" << eps_ << "}, {sigma:" << sigma_ << "}, {c6:" << c6_ << "}, {c12:" << c12_ << "}\n";
    }
    // X and Y are real positions, XS and YS are scaled positions
    virtual void CalcPotential(double* x,
                               double* xs,
                               double* y,
                               double* ys,
                               double* fpote) {
        std::fill(fpote, fpote + ndim_ + 1, 0.0);
     
        double rx[3], rxs[3];
        for (int i = 0; i < ndim_; ++i) {
            rx[i] = x[i] - y[i];
            rxs[i] = xs[i] - ys[i];
        }
        // Apply periodic boundary conditions
        periodic_boundary_conditions(space_->n_periodic, space_->unit_cell, space_->unit_cell_inv, rx, rxs);
        double rsq = 0.0;
        for (int i = 0; i < ndim_; ++i) {
            rsq += rx[i]*rx[i];
        }

        if (rsq < rcut2_) {
            double ffac, r6, rinv;

            rinv = 1.0/rsq;
            r6 = rinv*rinv*rinv;

            ffac = (12.0*c12_*r6 - 6.0*c6_)*r6*rinv;
            for (int i = 0; i < ndim_; ++i) {
                fpote[i] = rx[i]*ffac;
            }
            fpote[ndim_] = r6*(c12_*r6 - c6_);
        }
    }
};

#endif

