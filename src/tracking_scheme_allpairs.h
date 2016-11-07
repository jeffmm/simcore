#ifndef _SIMCORE_TRACKING_SCHEME_ALLPAIRS_H_
#define _SIMCORE_TRACKING_SCHEME_ALLPAIRS_H_

#include "tracking_scheme.h"

class TrackingSchemeAllPairs : public TrackingScheme {

  public:
    
    TrackingSchemeAllPairs() {
      name_ = "AllPairs";
    }
    virtual ~TrackingSchemeAllPairs() {}

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
    void GenerateAllPairs();
    void GenerateAllPairsSymmetric();
    void GenerateAllPairsAsymmetric();
};

#endif
