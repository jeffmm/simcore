#ifndef _SIMCORE_UBERENGINE_H_
#define _SIMCORE_UBERENGINE_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "interaction.h"
#include "interaction_engine.h"
//#include "kmc_engine.h"
#include "particle_engine.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <chrono>

class UberEngine {
  private:
    int ndim_,
        nperiodic_,
        nthreads_;
    system_parameters *params_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;
    al_set *anchors_;
    rng_properties rng_;

    ParticleEngine pengine_;
    InteractionEngine fengine_;
    //kmcEngine kengine_;

    std::vector<interaction_t> interactions_;

    // Statistics and timing
    std::chrono::time_point<std::chrono::high_resolution_clock> last_time_;
    std::chrono::time_point<std::chrono::high_resolution_clock> this_time_;
    int ndatapoints_;
    int ninteractions_;


  public:
    UberEngine() {}
    ~UberEngine() {
      PrintStatistics();
    }

  public:
    void DumpAll();
    void GenerateStatistics(int istep);
    void Init(system_parameters *pParams,
              space_struct *pSpace,
              std::vector<SpeciesBase*> *pSpecies,
              al_set *pAnchors,
              long seed);
    void Interact();
    void PrepOutputs();
    void PrintStatistics();
    //void StepKMC();
    void WriteOutputs(int istep);
};

#endif // _SIMCORE_UBERENGINE_H_
