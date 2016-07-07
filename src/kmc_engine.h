#ifndef _SIMCORE_KMC_ENGINE_H_
#define _SIMCORE_KMC_ENGINE_H_

#include "auxiliary.h"
#include "particle_tracking.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class kmcEngine {
  public:
    kmcEngine() {}
    ~kmcEngine() {}

    void Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, ParticleTracking *pTracking);
    void StepKMC();
    void PrepKMC();

  private:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    
    space_struct *space_;
    std::vector<SpeciesBase*>* species_;
};

#endif
