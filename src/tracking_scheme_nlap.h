#ifndef _SIMCORE_TRACKING_SCHEME_NEIGHBOR_LIST_ALLPAIRS_H_
#define _SIMCORE_TRACKING_SCHEME_NEIGHBOR_LIST_ALLPAIRS_H_

#include "tracking_scheme.h"

class TrackingSchemeNeighborListAllPairs : public TrackingScheme {
  
  public:

    TrackingSchemeNeighborListAllPairs() {
      name_ = "NeighborListAllPairs";
    }
    virtual ~TrackingSchemeNeighborListAllPairs() {}

    virtual void GenerateInteractions(bool pForceUpdate = false);
    virtual void Init(int pModuleID,
                      space_struct *pSpace,
                      PotentialBase *pPotentialBase,
                      std::vector<interaction_t> *pInteractions,
                      std::vector<SpeciesBase*> *pSpecies,
                      std::vector<Simple*> *pSimples,
                      std::unordered_map<int, int> *pOIDMap,
                      YAML::Node *pNode);
    virtual void Print();
    virtual void PrintStatistics();

  protected:

    virtual void CreateTrackingScheme();
    virtual void LoadSimples();
    void UpdateNeighborList();

    bool nl_update_ = false;
    double rcut_ = 0.0;
    double rcut2_ = 0.0;
    double skin_ = 0.0;
    double skin2_ = 0.0;
    double rcs2_ = 0.0;
    double half_skin2_ = 0.0;
};

#endif
