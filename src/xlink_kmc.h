#ifndef _SIMCORE_XLINK_KMC_H_
#define _SIMCORE_XLINK_KMC_H_

#include "auxiliary.h"
#include "kmc_base.h"

class Xlink;

class XlinkKMC : public KMCBase {
  protected:
    double eps_eff_0_1_[2];
    double eps_eff_1_2_[2];
    double on_rate_0_1_[2];
    double on_rate_1_2_[2];
    double alpha_;
    double barrier_weight_;
    double k_stretch_;
    double mrcut_;
    double mrcut2_;
    double velocity_;
    double max_length_;
    double rcutoff_1_2_;
    double r_equil_;

    int nsimples_;

    std::ofstream kmc_file;

    std::vector<Simple*>* simples_;
    std::unordered_map<int, int>* oid_position_map_;
    nl_list *neighbors_;

    void KMC_0_1();
    void KMC_1_0();
    void Update_0_1(Xlink *xit);
    void Update_1_2(Xlink *xit);
    void UpdateStage1(Xlink *xit, int *nbound1);

    void CalcCutoff();

  public:
    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking,
        SpeciesBase *spec1, SpeciesBase *spec2, int ikmc, YAML::Node &node,
        long seed);
    virtual void Print();
    virtual void Dump();
    virtual double GetMaxRcut() {
      return rcutoff_1_2_;
    }

    virtual void PrepKMC();
    virtual void StepKMC();
    virtual void UpdateKMC();

    virtual void PrepOutputs();
    virtual void WriteOutputs();

};

#endif
