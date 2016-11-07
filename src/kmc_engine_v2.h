#ifndef _SIMCORE_KMC_ENGINE_V2_H_
#define _SIMCORE_KMC_ENGINE_V2_H_

#include "auxiliary.h"
#include "helpers.h"
#include "interaction.h"
#include "kmc_base_v2.h"
#include "particle_engine.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class kmcEngineV2 {
  public:
    kmcEngineV2() {}
    ~kmcEngineV2() {}

    void Dump();
    void PreGenerateNeighbors();
    void Init(system_parameters *pParams,
              space_struct *pSpace,
              ParticleEngine *pTrackEngine,
              std::vector<interaction_t> *pInteractions,
              long seed);
    void InitMP();
    void PrepOutputs();
    void Print();
    void StepKMC();
    void WriteOutputs(int istep);

  protected:
      int ndim_;
      int nperiodic_;
      int nthreads_;
      int nsimples_;
      int nspecies_;
      int nkmcs_;

      system_parameters *params_;
      space_struct *space_;
      ParticleEngine *ptrack_;
      std::vector<interaction_t> *interactions_;
      std::vector<SpeciesBase*> *species_;
      std::vector<Simple*> *simples_;
      rng_properties rng_;

      rfh::factory kmc_factory_;
      std::vector<KMCBaseV2*> kmc_modules_;

      void AttachParticleEngine();
      void ParseKMC();
      void RegisterKMC();
};

#endif
