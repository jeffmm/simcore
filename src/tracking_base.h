// Force substructure base
// IE Cell lists, neighbor lists, neighbor lists with cell lists

#ifndef _SIMCORE_TRACKING_BASE_H_
#define _SIMCORE_TRACKING_BASE_H_

#include "auxiliary.h"
#include "neighbor_list_generic.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class TrackingBase {
  public:

    TrackingBase() {}
    virtual ~TrackingBase() {}

    virtual void Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, std::vector<Simple*> *pSimples, double pSkin);

    virtual void print() = 0; // print information
    virtual void dump() = 0;
    virtual void CreateSubstructure(double pRcut, nl_list** pNeighbors) = 0; // Create the substructure - depends on type
    virtual void UpdateTracking(bool pForceUpdate = false) = 0; // Update the tracking information for whatever we are

    std::string Name() const {return name_;}
    const int NUpdates() {return nupdates_;}

  protected:

    // Inputs that are needed
    int ndim_;
    int nperiodic_;
    int nupdates_;
    double skin_ = 0.0;
    double rcut_ = 0.0;
    double box_[3];

    space_struct* space_;
    std::vector<Simple*> *simples_;
    std::vector<SpeciesBase*> *species_;
    nl_list* neighbors_;

    // Derived quantities
    int nthreads_;
    int nsimples_;
    std::string name_ = "TrackingBase";

};

template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<TrackingBase, T>::value>::type>
T* trackingFactory(ARGS&&... args) {
    T* mtrack{ new T{std::forward<ARGS>(args)...} };
    
    return mtrack;
}

#endif
