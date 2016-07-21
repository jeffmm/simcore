#ifndef _SIMCORE_UBERENGINE_H_
#define _SIMCORE_UBERENGINE_H_

#include "auxiliary.h"
#include "interaction_engine.h"
#include "kmc_engine.h"
#include "minimum_distance.h"
#include "potential_manager.h"
#include "particle_tracking.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class UberEngine {
  private:
    int n_dim_,
        n_periodic_,
        max_overlap_,
        draw_flag_,
        nthreads_;
    bool draw_;
    double dr_[3],
           contact1_[3],
           contact2_[3],
           dr_mag_,
           dr_mag2_,
           buffer_mag_,
           buffer_mag2_,
           skin_;
    system_parameters *params_;
    FTYPE force_type_;
    space_struct *space_;
    std::vector<SpeciesBase*> *species_;
    PotentialManager potentials_;
    ParticleTracking ptrack_;
    InteractionEngine fengine_; //fengine = force engine.  get it?
    kmcEngine kengine_;
    rng_properties rng_;

  public:
    UberEngine() {}
    ~UberEngine() {}

    std::vector<graph_struct> draw_array_;
  public:
    void Init(system_parameters *pParams, space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, long seed);
    void InteractMP();
    void StepKMC();
    void DumpAll();
    void InitPotentials();
    void Draw(std::vector<graph_struct*> * graph_array);
    void PrepOutputs();
    void WriteOutputs();
};

#endif // _SIMCORE_UBERENGINE_H_
