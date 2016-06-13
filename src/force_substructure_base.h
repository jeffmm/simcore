// Force substructure base
// IE Cell lists, neighbor lists, neighbor lists with cell lists

#ifndef _SIMCORE_FORCE_SUBSTRUCTURE_BASE_H_
#define _SIMCORE_FORCE_SUBSTRUCTURE_BASE_H_

#include "auxiliary.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class ForceSubstructureBase {
  public:

    ForceSubstructureBase() {}
    virtual ~ForceSubstructureBase() {}

    virtual void Init(space_struct *pSpace, double pSkin);
    virtual void LoadFlatSimples(std::vector<Simple*> pSimples);

    virtual void CreateSubstructure(double pRcut) = 0; // Create the substructure - depends on type

  protected:

    // Inputs that are needed
    int ndim_;
    int nperiodic_;
    double skin_;
    double box_[3];

    space_struct* space_;
    std::vector<Simple*> simples_;

    // Derived quantities
    int nthreads_;
    int nparticles_;

};

#endif
