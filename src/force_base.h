// Basic type of force, will depend on what is the underlying force generation structure
// i.e. brute force, cell list, neighbor list, etc

#ifndef _CYTOSCORE_FORCE_BASE_H_
#define _CYTOSCORE_FORCE_BASE_H_

#include "auxiliary.h"
#include "species.h"
#include "minimum_distance.h"
#include "potential_manager.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class ForceBase {
  public:

    ForceBase() {}
    virtual ~ForceBase() {}

    virtual void Init(space_struct *pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);
    virtual void InitPotentials(std::vector<SpeciesBase*> pSpecies);
    virtual void MinimumDistance(Simple* o1, Simple* o2, interactionmindist& idm);

    virtual void InitMP() = 0; // init the underlying scheme
    virtual void UpdateScheme() = 0; // updates the underlying scheme (cell, neighbor, etc)
    virtual void Interact() = 0; // Main interaction routine!!!!!!!
    
  protected:

    // Inputs that we need
    int nsys_;
    int ndim_;
    int nperiodic_;
    double skin_;
    double box_[3];

    // Derived/calculated quantities
    int nthreads_;
    int nparticles_;
    
    space_struct* space_;

    std::vector<Simple*> simples_;
    std::vector<double> frc_; // Force superarray for threading
    std::vector<double> prc_energy_; // Energy superarray for threading

    // Classes for managers
    PotentialManager potentials_;
};


// Force factory using templates
template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<ForceBase, T>::value>::type>
T* forceFactory(ARGS&&... args) {
    T* mforce{ new T{std::forward<ARGS>(args)...} };
    
    return mforce;
}

#endif
