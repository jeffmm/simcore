#ifndef _SIMCORE_KMC_BASE_H_
#define _SIMCORE_KMC_BASE_H_

#include "auxiliary.h"
#include "particle_tracking.h"
#include "potential_manager.h"
#include "species.h"

// KMC base for all the kmc routines that we can use
class KMCBase {
  protected:
    int ndim_;
    int nperiodic_;

    rng_properties rng_;
    space_struct *space_;
    ParticleTracking *tracking_;
    PotentialManager *potentials_;
    SID sid1_;
    SID sid2_;
    SpeciesBase *spec1_;
    SpeciesBase *spec2_;

  public:
    KMCBase() {}
    virtual ~KMCBase() {}

    virtual void Init(space_struct *pSpace,
                      ParticleTracking *pTracking,
                      PotentialManager *pPotentials,
                      SpeciesBase *spec1,
                      SpeciesBase *spec2,
                      int ikmc, YAML::Node &node, long seed) {
      space_ = pSpace;
      ndim_ = space_->n_dim;
      nperiodic_ = space_->n_periodic;
      rng_.init(seed);
      tracking_ = pTracking;
      std::string sid1s   = node["kmc"][ikmc]["sid1"].as<std::string>();
      std::string sid2s   = node["kmc"][ikmc]["sid2"].as<std::string>();
      sid1_ = StringToSID(sid1s);
      sid2_ = StringToSID(sid2s);
      spec1_ = spec1;
      spec2_ = spec2;
      potentials_ = pPotentials;
    }
    virtual double GetMaxRcut() {return 0.0;}
    virtual void PrepKMC() {}
    virtual void StepKMC() {}
    virtual void UpdateKMC() {}
    virtual void TransferForces() {}

    virtual void PrepOutputs() {}
    virtual void WriteOutputs(int istep) {}

    virtual void Print() {
      std::cout << "\tsids: [" << (int)sid1_ << ", " << (int)sid2_ << "]\n";
    }
    virtual std::pair<SID, SID> GetSIDs() {return std::make_pair(sid1_, sid2_);}
    virtual void Dump() {
      printf("ERROR, need to override this function!\n");
      exit(1);
    }

    void SetRNGState(const std::string& filename) {
      // Load the rng state from binary file
      FILE* pfile = fopen(filename.c_str(), "r");
      auto retval = gsl_rng_fread(pfile, rng_.r);
      if (retval != 0) {
        std::cout << "Reading rng state failed " << retval << std::endl;
      }
    }
};

typedef std::map<sid_pair, KMCBase*> kmc_map;

#endif // _SIMCORE_KMC_BASE_H_
