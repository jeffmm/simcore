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
    ~ParticleEngine();

    void Init(system_parameters *pParams,
              space_struct *pSpace,
              std::vector<SpeciesBase*> *pSpecies,
              al_set *pAnchors,
              std::vector<interaction_t> *pInteractions,
              long seed);

    void CheckTriggerUpdate();
    void CreateExternalPotential(YAML::Node *subnode, int potidx);
    TrackingScheme *CreateKMCTracking(YAML::Node *subnode);
    void CreateTracking();
    void Dump();
    void DumpInteractions();
    void DumpSimples();
    void LoadSimples();
    void Print();
    void PrintStatistics();
    void RegisterSchemes();
    void UpdateInteractions();

    std::vector<Simple*> *GetSimples() {return &simples_;}
    std::vector<SpeciesBase*> *GetSpecies() {return species_;}
    std::unordered_map<int, int> *GetOIDPositionMap() {return &oid_position_map_;}
    std::vector<interaction_t> *GetInteractions() {return interactions_;}

  protected:
    bool trigger_update_ = false;

    int ndim_;
    int nperiodic_;
    int nthreads_;
    int npots_;
    int nsys_;
    int nsimples_;
    
    system_parameters *params_;
    space_struct *space_;
    al_set* anchors_;
    rng_properties rng_;

    std::vector<interaction_t> *interactions_;
    std::vector<SpeciesBase*> *species_;
    std::vector<Simple*> simples_;
    std::unordered_map<int, int> oid_position_map_;

    // Things we actually own
    YAML::Node node_;
    PotentialManagerV2 potentials_;
    rfh::factory scheme_factory_;

    // Potential mappings
    std::vector<TrackingScheme*> tracking_;

};

#endif /* _SIMCORE_PARTICLE_ENGINE_H_ */
