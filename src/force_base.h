// Basic type of force, will depend on what is the underlying force generation structure
// i.e. brute force, cell list, neighbor list, etc

#ifndef _CYTOSCORE_FORCE_BASE_H_
#define _CYTOSCORE_FORCE_BASE_H_

#include <unordered_map>
#include <memory>

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
    virtual ~ForceBase() {
        space_ = NULL;
        simples_.clear();
        //frc_.clear();
        //prc_energy_.clear();
        delete[] frc_;
        delete[] trqc_;
        delete[] prc_energy_;
        oid_position_map_.clear();
    }

    virtual void Init(space_struct *pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);
    virtual void InitPotentials(std::vector<SpeciesBase*> pSpecies);
    virtual void MinimumDistance(Simple* o1, Simple* o2, interactionmindist& idm);

    virtual void Finalize() = 0; // AFter we do everything, we must finalize it for safety
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
    bool initialized_ = false;
    int nthreads_;
    int nparticles_;
    
    space_struct* space_;

    std::vector<Simple*> simples_;
    double* frc_; // Force superarray
    double* trqc_; // torque superarray
    double* prc_energy_; // Potential energy superarray
    std::unordered_map<int, int> oid_position_map_; // oid to position mapping!!!

    // Classes for managers
    PotentialManager potentials_;
};


// Force factory using templates
/*template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<ForceBase, T>::value>::type>
std::shared_ptr<T> forceFactory(ARGS&&... args) {
    std::shared_ptr<T> mforce{ new T{std::forward<ARGS>(args)...} };
    
    return mforce;
}*/
template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<ForceBase, T>::value>::type>
T* forceFactory(ARGS&&... args) {
    T* mforce{ new T{std::forward<ARGS>(args)...} };
    
    return mforce;
}

#endif
