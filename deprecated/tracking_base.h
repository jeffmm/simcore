// Force substructure base
// IE Cell lists, neighbor lists, neighbor lists with cell lists

#ifndef _SIMCORE_TRACKING_BASE_H_
#define _SIMCORE_TRACKING_BASE_H_

#include "auxiliary.h"
#include "neighbor_list_generic.h"
#include "potential_manager.h"
#include "species.h"

#include <unordered_set>

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class TrackingBase {
  public:

    TrackingBase() {}
    virtual ~TrackingBase() {
      if (unique_rids_)
        delete unique_rids_;
      if (rid_self_check_)
        delete rid_self_check_;
      if (rid_check_local_)
        delete[] rid_check_local_;
      //delete[] rid_interactions_;
    }

    virtual void Init(space_struct *pSpace,
                      std::vector<SpeciesBase*> *pSpecies,
                      std::vector<Simple*> *pSimples,
                      PotentialManager *pPotentials,
                      double pSkin);

    virtual void print() = 0; // print information
    virtual void dump() = 0;
    virtual void CreateSubstructure(double pRcut, nl_list** pNeighbors) = 0; // Create the substructure - depends on type
    virtual void UpdateTracking(bool pForceUpdate = false) = 0; // Update the tracking information for whatever we are
    virtual void Rebuild(nl_list **pNeighbors);
    virtual void UpdateRcut(double pRcut) = 0; // Force the change in rcut

    const std::string Name() const {
      std::cout << "My name is " << name_ << std::endl;
      return name_;
    }
    const int NUpdates() {return nupdates_;}

  protected:

    // Inputs that are needed
    int ndim_;
    int nperiodic_;
    int nupdates_;
    int nthreads_;
    int nsimples_;
    int nrigids_;
    int maxrigid_;
    double skin_ = 0.0;
    double rcut_ = 0.0;
    double box_[3];

    std::string name_ = "TrackingBase";

    space_struct* space_;
    std::vector<Simple*> *simples_;
    std::vector<SpeciesBase*> *species_;
    nl_list* neighbors_;
    PotentialManager *potentials_;
  
    // Rigid stuff
    std::unordered_set<std::pair<int, int>, hashh::pair_hash>* rid_interactions_;
    std::vector<bool>* rid_self_check_;
    //std::unordered_set<int>* unique_rids_;
    std::set<int>* unique_rids_;
    std::unordered_set<int>** rid_check_local_;


};

template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<TrackingBase, T>::value>::type>
T* trackingFactory(ARGS&&... args) {
    T* mtrack{ new T{std::forward<ARGS>(args)...} };
    
    return mtrack;
}

#endif
