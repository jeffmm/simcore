#ifndef _SIMCORE_PARTICLE_TRACKING_H_
#define _SIMCORE_PARTICLE_TRACKING_H_

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
    const int GetNSimples() {return nsimples_;}
    std::vector<Simple*>* GetSimples() {return &simples_;}

    void Print();
    void Dump();

  private:
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

    // Everybody uses a neighbor list now
    // It's just how we build it that's different
    nl_list* neighbors_;
    TrackingBase* tracking_;
};

#endif
