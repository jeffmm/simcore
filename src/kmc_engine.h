#ifndef _SIMCORE_KMC_ENGINE_H_
#define _SIMCORE_KMC_ENGINE_H_

#include "auxiliary.h"
#include "helpers.h"
#include "kmc_base.h"
#include "particle_tracking.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class kmcEngine {
  public:
    kmcEngine() {}
    ~kmcEngine() {}

    void Init(space_struct *pSpace,
              std::vector<SpeciesBase*> *pSpecies,
              ParticleTracking *pTracking,
              long seed,
              char *pFname);
    void InitMP();
    void RunKMC();
    void StepKMC();
    void PrepKMC();
    void UpdateKMC();
    void Print();
    void Dump();

    void AddKMC(SID sid1, SID sid2, KMCBase *kmc);
    void RegisterKMC();
    void ParseKMC();

  private:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int nsimples_;
    int nkmcs_;

    std::string fname_;
    
    space_struct *space_;
    std::vector<SpeciesBase*>* species_;
    std::vector<Simple*>* simples_;
    ParticleTracking *tracking_;
    rng_properties rng_;

    rfh::factory kmc_factory_;
    kmc_map kmc_map_;
};

#endif
