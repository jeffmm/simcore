// Base class for all potentials
// All potentials need to know about the cutoff, particle 1, particle 2
// any other information about interactions...
#ifndef BUFFMD_POTENTIAL_BASE_H_
#define BUFFMD_POTENTIAL_BASE_H_

#include <array>

#include "helpers.h"

class PotentialBase {
public:
    PotentialBase(double pRcut, double pBox) : rcut_(pRcut), box_(pBox) {
        rcsq_ = rcut_ * rcut_;
        boxby2_ = 0.5 * box_;
    }

    virtual ~PotentialBase() {

    }
    
    virtual void CalcPotential(double x[3],
                               double y[3],
                               double* fpote) {
        buffmd::unused(x,y);
        std::fill(fpote, fpote + 4, 0.0);
    }

    virtual void print() {
        printf("\t{rcut:%f}, {rcsq:%f}, {box:%f}\n", rcut_, rcsq_, box_);
    }

protected:
    double rcut_, rcsq_;
    double box_, boxby2_;
};

// Potential factory using templates, similar to particle factory
template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<PotentialBase, T>::value>::type>
T* potentialFactory(ARGS&&... args) {
    T* pot{ new T{ std::forward<ARGS>(args)...} };

    return pot;
}
#endif
