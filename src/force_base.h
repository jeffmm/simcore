// Basic type of force, will depend on what is the underlying force generation structure
// i.e. brute force, cell list, neighbor list, etc

#ifndef _SIMCORE_FORCE_BASE_H_
#define _SIMCORE_FORCE_BASE_H_

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
        delete[] frc_;
        delete[] trqc_;
        delete[] prc_energy_;
        oid_position_map_.clear();
    }

    virtual void Init(space_struct *pSpace, double pSkin);
    virtual void LoadSimples(std::vector<SpeciesBase*> pSpecies);
    virtual void InitPotentials(PotentialManager *pPotentials);

    // IO routines of awfulness
    virtual void print();
    virtual void printSpecifics() = 0; // must be overridden
    virtual void dump();

    void InteractParticlesMP(Simple *part1, Simple* part2, double **fr, double **tr, double *pr_energy);
    void ReduceParticlesMP();

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
    int nthreads_ = 1;
    int nparticles_ = 0;
    double max_rcut_ = 0.0;
    std::string name_ = "ForceBase";
    
    space_struct* space_;

    std::vector<Simple*> simples_;
    double* frc_; // Force superarray
    double* trqc_; // torque superarray
    double* prc_energy_; // Potential energy superarray
    std::unordered_map<int, int> oid_position_map_; // oid to position mapping!!!

    // Classes for managers
    PotentialManager *potentials_;
};


template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<ForceBase, T>::value>::type>
T* forceFactory(ARGS&&... args) {
    T* mforce{ new T{std::forward<ARGS>(args)...} };
    
    return mforce;
}

#endif
