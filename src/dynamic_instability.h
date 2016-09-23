#ifndef _SIMCORE_DYNAMIC_INSTABILITY_H_
#define _SIMCORE_DYNAMIC_INSTABILITY_H_

#include "auxiliary.h"
#include "br_rod.h"
#include "kmc_base.h"

class DynamicInstabilityKMC : public KMCBase {

  protected:
    double delta_;
    // 6 transitions
    double f_shrink_to_grow_;
    double f_shrink_to_pause_;
    double f_pause_to_shrink_;
    double f_pause_to_grow_;
    double f_grow_to_pause_;
    double f_grow_to_shrink_;
    // 6 transitions as probabilities
    double p_stg_,
           p_stp_,
           p_ptg_,
           p_pts_,
           p_gtp_,
           p_gts_;
    // velocities
    double v_poly_;
    double v_depoly_;
    // Max and min lengths
    double max_length_;
    double min_length_;

  public:

    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking,
        PotentialManager *pPotentials,
        SpeciesBase *spec1, SpeciesBase *spec2, int ikmc, YAML::Node &node,
        long seed);
    virtual void Print();
    virtual void Dump();
    virtual double GetMaxRcut() {return 0.0;}

    virtual void StepKMC();

    // Specific Functions
    void UpdatePolymerizationState();
    void GrowBonds();
};

#endif
