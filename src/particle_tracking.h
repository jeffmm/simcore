#ifndef _SIMCORE_PARTICLE_TRACKING_H_
#define _SIMCORE_PARTICLE_TRACKING_H_

#include <unordered_map>

#include "auxiliary.h"
#include "species.h"
#include "minimum_distance.h"
#include "neighbor_list_generic.h"
#include "potential_manager.h"
#include "tracking_base.h"
#include "tracking_allpairs.h"
#include "tracking_neighborlist_ap.h"

class ParticleTracking {
  public:
    ParticleTracking() {}
    ~ParticleTracking() {
      delete[] neighbors_;
      oid_position_map_.clear();
      printf("********\n");
      printf("ParticleTracking Stats ->\n");
      printf("\tUpdates: %d\n", tracking_->NUpdates());
    }

    void Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, double pSkin, FTYPE pFtype);
    void LoadSimples();
    void InitPotentials(PotentialManager *pPotentials);
    void InitTracking();
    void CheckOverlaps(int pMaxOverlaps);

    void UpdateTracking(bool pForceUpdate = false);

    // Getters and setters
    nl_list* GetNeighbors() {return neighbors_;}
    std::unordered_map<int, int>* GetOIDPositionMap() {return &oid_position_map_;}
    const int GetNSimples() {return nsimples_;}
    std::vector<Simple*>* GetSimples() {return &simples_;}
    const bool TriggerUpdate() {return trigger_update_;}

    void Print();
    void Dump();

  private:
    bool trigger_update_;
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int nsys_;
    int nsimples_;
    double rcut_;
    double skin_;
    double box_[3];
    FTYPE ftype_;
    
    space_struct *space_;
    PotentialManager *potentials_;
    std::vector<Simple*> simples_;
    std::vector<SpeciesBase*>* species_;
    std::unordered_map<int, int> oid_position_map_; // oid to position mapping!!!

    // Everybody uses a neighbor list now
    // It's just how we build it that's different
    nl_list* neighbors_;
    TrackingBase* tracking_;
};

#endif
