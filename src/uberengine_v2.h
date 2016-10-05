#ifndef _SIMCORE_UBERENGINE_V2_H_
#define _SIMCORE_UBERENGINE_V2_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "interaction.h"
#include "interaction_engine_v2.h"
#include "particle_engine.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class UberEngineV2 {
  private:
    int ndim_,
        nperiodic_,
        nthreads_;
    system_parameters *params_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;
    al_set *anchors_;
    rng_properties rng_;

    // v2.0 stuff
    ParticleEngine pengine_;
    InteractionEngineV2 fengine_;

    std::vector<interaction_t> interactions_;

  public:
    UberEngineV2() {}
    ~UberEngineV2() {}

  public:
    void DumpAll();
    void Init(system_parameters *pParams,
              space_struct *pSpace,
              std::vector<SpeciesBase*> *pSpecies,
              al_set *pAnchors,
              long seed);
    void InteractMP();
};

#endif // _SIMCORE_UBERENGINE_V2_H_
