#ifndef _CYTOSCORE_LJ126_H_
#define _CYTOSCORE_LJ126_H_

#include <cmath>

#include "potential_base.h"

class LJ126 : public PotentialBase {
  protected:
    double eps_, sigma_;
    double c12_, c6_;
  public:
    LJ126(int pNdim, double pEps, double pSigma, double pRcut, double pBox) : PotentialBase(pNdim, pRcut, pBox), eps_(pEps), sigma_(pSigma) {
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }
    virtual void print() {
        std::cout << "Lennard-Jones 12-6 potential:\n";
        PotentialBase::print();
        std::cout << "\t{eps:" << eps_ << "}, {sigma:" << sigma_ << "}, {c6:" << c6_ << "}, {c12:" << c12_ << "}\n";
    }
    virtual void CalcPotential(double* x,
                               double* y,
                               double* fpote) {
        std::fill(fpote, fpote + ndim_ + 1, 0.0);
        
        //double rx = buffmd::pbc(x[0] - y[0], boxby2_, box_);
        //double ry = buffmd::pbc(x[1] - y[1], boxby2_, box_);
        //double rz = buffmd::pbc(x[2] - y[2], boxby2_, box_);
        double rsq = rx*rx + ry*ry + rz*rz;

        if (rsq < rcut2_) {
            double ffac, r6, rinv;

            rinv = 1.0/rsq;
            r6 = rinv*rinv*rinv;

            ffac = (12.0*c12_*r6 - 6.0*c6_)*r6*rinv;
            fpote[0] = rx*ffac;
            fpote[1] = ry*ffac;
            fpote[2] = rz*ffac;
            fpote[3] = r6*(c12_*r6 - c6_);
        }
    }
};

#endif

