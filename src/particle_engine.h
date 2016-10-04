#ifndef _SIMCORE_PARTICLE_ENGINE_H_
#define _SIMCORE_PARTICLE_ENGINE_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "interaction.h"
#include "potential_manager_v2.h"
#include "tracking_scheme.h"
#include "species.h"

#include <yaml-cpp/yaml.h>

class ParticleEngine {
  
  public:

    ParticleEngine() {}
    ~ParticleEngine() {}

    void Init(system_parameters *pParams,
              space_struct *pSpace,
              std::vector<SpeciesBase*> *pSpecies,
              al_set *pAnchors,
              std::vector<interaction_t> *pInteractions,
              long seed);

    void CreateTracking();
    void CreateExternalPotential(YAML::Node *subnode, int potidx);
    void Print();
    void RegisterSchemes();

  private:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int npots_;
    
    system_parameters *params_;
    space_struct *space_;
    al_set* anchors_;
    rng_properties rng_;

    std::vector<interaction_t> *interactions_;
    std::vector<SpeciesBase*> *species_;

    // Things we actually own
    YAML::Node node_;
    PotentialManagerV2 potentials_;
    rfh::factory scheme_factory_;

    // Potential mappings
    std::vector<TrackingScheme*> tracking_;
};

#endif /* _SIMCORE_PARTICLE_ENGINE_H_ */
