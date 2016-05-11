// Lennard Jones 12-6 potential

#ifndef BUFFMD_LJ126_H_
#define BUFFMD_LJ126_H_

#include <cmath>

#include "potential_base.h"
#include "helpers.h"

class LJ126 : public PotentialBase {
public:
    LJ126(double pEps, double pSigma, double pRcut, double pBox) : PotentialBase(pRcut, pBox), eps_(pEps), sigma_(pSigma) {
        c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
        c6_  = 4.0 * eps_ * pow(sigma_,  6.0);
    }

    virtual void print() {
        printf("Lennard-Jones 12-6 potential:\n");
        PotentialBase::print();
        printf("\t{eps:%f}, {sigma:%f}, {c6:%f}, {c12:%f}\n", eps_, sigma_, c6_, c12_);
    }

    virtual void CalcPotential(double x[3],
                               double y[3],
                               double* fpote) {
        std::fill(fpote, fpote + 4, 0.0);
        
        double rx = buffmd::pbc(x[0] - y[0], boxby2_, box_);
        double ry = buffmd::pbc(x[1] - y[1], boxby2_, box_);
        double rz = buffmd::pbc(x[2] - y[2], boxby2_, box_);
        double rsq = rx*rx + ry*ry + rz*rz;

        if (rsq < rcsq_) {
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

protected:
    double eps_, sigma_;
    double c12_, c6_;
};

#endif

