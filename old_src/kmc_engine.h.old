#ifndef _SIMCORE_KMC_ENGINE_H_
#define _SIMCORE_KMC_ENGINE_H_

#include "auxiliary.h"
#include "helpers.h"
#include "kmc_base.h"
#include "particle_tracking.h"
#include "potential_manager.h"
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
    void InitPotentials(PotentialManager *pPotentials);
    void InitMP();
    void RunKMC();
    void Print();
    void Dump();
    void PrepOutputs();
    void WriteOutputs(int istep);

    void RegisterKMC();
    void ParseKMC();

    double GetMaxRcut();

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
    PotentialManager *potentials_;
    rng_properties rng_;

    rfh::factory kmc_factory_;
    std::vector<KMCBase*> kmc_modules_;
};

#endif
