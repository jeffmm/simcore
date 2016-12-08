#ifndef _SIMCORE_PARTICLE_ENGINE_H_
#define _SIMCORE_PARTICLE_ENGINE_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "interaction.h"
#include "potential_manager.h"
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
    TrackingScheme *CreateKMCExternal(YAML::Node *subnode);
    PotentialBase *CreateKMCInternal(YAML::Node *subnode);
    void CreateInternalPotential(YAML::Node *subnode, int potidx);
    void CreateBoundaryPotential(YAML::Node *subnode, int potidx);
    void CreateTetherPotential(YAML::Node *subnode, int potidx);
    void CreateTracking();
    void Dump();
    void DumpInteractions();
    void DumpNeighbors();
    void DumpSimples();
    void LoadSimples();
    void Print();
    void PrintPotentials();
    void PrintStatistics();
    void RegisterSchemes();
    void UpdateInteractions(bool pForceUpdate = false);

    std::vector<Simple*> *GetSimples() {return &simples_;}
    std::vector<SpeciesBase*> *GetSpecies() {return species_;}
    std::unordered_map<int, int> *GetOIDPositionMap() {return &oid_position_map_;}
    std::vector<interaction_t> *GetInteractions() {return interactions_;}
    al_set *GetAnchors() {return anchors_;}
    const bool GetTriggerUpdate() {return trigger_update_;}

    void SetTriggerUpdate(bool trigger_update) {trigger_update_=trigger_update;}

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
    std::vector<interaction_t> m_internal_interactions_; // internal interactions that never change
    std::vector<interaction_t> m_boundary_interactions_; // boundary interactions that never change
    std::vector<interaction_t> m_tether_interactions_; // tether interactions that never change
    std::vector<SpeciesBase*> *species_;
    std::vector<Simple*> simples_;
    std::unordered_map<int, int> oid_position_map_;

    // Things we actually own
    YAML::Node node_;
    PotentialManager potentials_;
    rfh::factory scheme_factory_;

    // Potential mappings
    std::vector<TrackingScheme*> tracking_;

};

#endif /* _SIMCORE_PARTICLE_ENGINE_H_ */
