#ifndef _SIMCORE_BR_BINDUNBIND_H_
#define _SIMCORE_BR_BINDUNBIND_H_

#include "auxiliary.h"
#include "kmc_base.h"

class BrBindUnbind : public KMCBase {
  protected:
    double eps_eff_;
    double on_rate_;
    double alpha_;
    double mrcut_;

    int nsimples_;
    std::vector<Simple*>* simples_;

    void Bind();
    void Unbind();

  public:
    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking,
        SpeciesBase *spec1, SpeciesBase *spec2, int ikmc, YAML::Node &node,
        long seed);
    virtual void Print();
    virtual void Dump();

    virtual void PrepKMC();
    virtual void StepKMC();
    virtual void UpdateKMC();
};

#endif
