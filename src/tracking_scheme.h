#ifndef _SIMCORE_TRACKING_SCHEME_H_
#define _SIMCORE_TRACKING_SCHEME_H_

#include "auxiliary.h"
#include "interaction.h"
#include "neighbor_list.h"
#include "potential_base.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <unordered_set>
#include <set>

class TrackingScheme {
  public:

    TrackingScheme() {}
    virtual ~TrackingScheme() {
      if (unique_rids_)
        delete unique_rids_;
      if (rid_self_check_)
        delete rid_self_check_;
      if (rid_check_local_)
        delete[] rid_check_local_;
      if (mneighbors_)
        delete[] mneighbors_;
    }

    // Virtual functions
    virtual void Dump();
    virtual void GenerateInteractions(bool pForceUpdate = false) = 0;
    virtual void GenerateStatistics();
    virtual void Init(int pModuleID,
                      space_struct *pSpace,
                      PotentialBase *pPotentialBase,
                      std::vector<interaction_t> *pInteractions,
                      std::vector<SpeciesBase*> *pSpecies,
                      std::vector<Simple*> *pSimples,
                      std::unordered_map<int, int>* pOIDMap,
                      YAML::Node *pNode);
    virtual void Print();
    virtual void PrintStatistics() = 0;

    // Non virtual functions
    ptype GetModuleType() {return type_;}
    int GetModuleID() {return moduleid_;}

  protected:

    int ndim_ = -1;
    int nperiodic_ = -1;
    int nthreads_ = -1;
    int nsimples_ = -1;
    int nmsimples_ = -1;
    int nupdates_ = 0;
    int maxrigid_ = 0;
    int moduleid_ = -1;

    SID sid0_;
    SID sid1_;
    SID kmc_target_;
    ptype type_;
    SpeciesBase *spec0_ = nullptr;
    SpeciesBase *spec1_ = nullptr;

    std::string name_ = "TrackingScheme";

    space_struct *space_;
    PotentialBase *pbase_;

    std::vector<interaction_t> *interactions_;
    std::vector<interaction_t> m_interactions_; // internal interactions if no update needed
    std::vector<Simple*> *simples_;
    std::vector<Simple*> m_simples_;
    std::unordered_map<int, int> *oid_position_map_;

    void CreateKMCNeighbors();
    virtual void LoadSimples();

    virtual void CreateTrackingScheme() = 0;

    // Rigid stuff
    std::vector<bool> *rid_self_check_;
    std::set<int> *unique_rids_;
    std::unordered_set<int>** rid_check_local_;

    // Statistics and timing
    std::chrono::time_point<std::chrono::high_resolution_clock> last_time_;
    std::chrono::time_point<std::chrono::high_resolution_clock> this_time_;
    double avg_update_time_;
    double avg_occupancy_;

    // KMC stuff
    nl_kmc_list* mneighbors_ = nullptr;

};

#endif
