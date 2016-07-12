#ifndef _SIMCORE_KMC_BASE_H_
#define _SIMCORE_KMC_BASE_H_

#include "auxiliary.h"
#include "particle_tracking.h"
#include "species.h"

// KMC base for all the kmc routines that we can use
class KMCBase {
  protected:
    int ndim_;

    rng_properties rng_;
    space_struct *space_;
    ParticleTracking *tracking_;
    SID sid1_;
    SID sid2_;

  public:
    KMCBase() {}
    virtual ~KMCBase() {}

    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking, int ikmc, YAML::Node &node, long seed) {
      space_ = pSpace;
      ndim_ = space_->n_dim;
      rng_.init(seed);
      tracking_ = pTracking;
      std::string sid1s   = node["kmc"][ikmc]["sid1"].as<std::string>();
      std::string sid2s   = node["kmc"][ikmc]["sid2"].as<std::string>();
      sid1_ = StringToSID(sid1s);
      sid2_ = StringToSID(sid2s);
    }
    virtual void RunKMC(SpeciesBase *spec1, SpeciesBase *spec2) {}
    virtual void Print() {
      printf("\tsids: [%d, %d]\n", sid1_, sid2_);
    }
};

typedef std::map<sid_pair, KMCBase*> kmc_map;

#endif // _SIMCORE_KMC_BASE_H_
