#ifndef _SIMCORE_KMC_BASE_V2_H_
#define _SIMCORE_KMC_BASE_V2_H_

#include "auxiliary.h"
#include "interaction.h"
#include "neighbor_list.h"
#include "particle_engine.h"
#include "species.h"

// KMC base for all the kmc routines that we can use
class KMCBaseV2 {
  protected:
    int ndim_;
    int nperiodic_;

    rng_properties rng_;
    space_struct *space_;
    ParticleEngine *ptrack_;
    SID sid1_;
    SID sid2_;
    SpeciesBase *spec1_;
    SpeciesBase *spec2_;

    std::vector<std::vector<neighbor_kmc_t>> nl_kmc_;
    std::vector<Simple*> *simples_;
    std::vector<SpeciesBase*> *species_;
    std::vector<interaction_t> *interactions_;
    std::unordered_map<int, int> *oid_position_map_;

  public:
    KMCBaseV2() {}
    virtual ~KMCBaseV2() {}

    virtual void Init(space_struct *pSpace,
                      ParticleEngine *pTrackEngine,
                      SpeciesBase *spec1,
                      SpeciesBase *spec2,
                      YAML::Node *subnode,
                      long seed) {
      YAML::Node node = *subnode;
      space_ = pSpace;
      ndim_ = space_->n_dim;
      nperiodic_ = space_->n_periodic;
      rng_.init(seed);
      ptrack_ = pTrackEngine;
      std::string sid1s = node["sid1"].as<std::string>();
      std::string sid2s = node["sid2"].as<std::string>();
      sid1_ = StringToSID(sid1s);
      sid2_ = StringToSID(sid2s);
      spec1_ = spec1;
      spec2_ = spec2;

      // Attach to the particle tracking
      simples_ = ptrack_->GetSimples();
      species_ = ptrack_->GetSpecies();
      oid_position_map_ = ptrack_->GetOIDPositionMap();
      interactions_ = ptrack_->GetInteractions();
    }
    virtual double GetMaxRcut() {return 0.0;}
    virtual void GenerateKMCNeighborList() {}
    virtual void PrepKMC() {}
    virtual void StepKMC() {}
    virtual void UpdateKMC() {}
    virtual void TransferForces() {}
    virtual void GenerateTrackingScheme(YAML::Node *subnode) {
      std::cout << "Should be overridden\n";
    }

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

#endif // _SIMCORE_KMC_BASE_V2_H_
