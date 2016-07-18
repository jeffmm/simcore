#ifndef _SIMCORE_BR_BINDUNBIND_H_
#define _SIMCORE_BR_BINDUNBIND_H_

#include "auxiliary.h"
#include "kmc_base.h"

class BrBindUnbind : public KMCBase {
  protected:
    double eps_eff_;
    double on_rate_;

    int nsimples_;
    std::vector<Simple*>* simples_;

    void Bind();
    void Unbind(SpeciesBase* spec);
    void FinishKMC(SpeciesBase* spec);

  public:
    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking, int ikmc, YAML::Node &node, long seed);
    virtual void RunKMC(SpeciesBase *spec1, SpeciesBase *spec2);
    virtual void Print();
};

#endif
