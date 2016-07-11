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

    void Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, ParticleTracking *pTracking, long seed);
    void InitMP();
    void StepKMC();
    void PrepKMC();
    void Dump();

    // Hardcoded function for now
    void HardcodedBindUnbind();
    void HardcodedBind();
    void HardcodedUnbind();

  private:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int nsimples_;
    
    space_struct *space_;
    std::vector<SpeciesBase*>* species_;
    std::vector<Simple*>* simples_;
    ParticleTracking *tracking_;
    rng_properties rng_;
};

#endif
