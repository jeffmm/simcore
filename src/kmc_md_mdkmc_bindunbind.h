#ifndef _SIMCORE_MD_MDKMC_BINDUNBIND_H_
#define _SIMCORE_MD_MDKMC_BINDUNBIND_H_

#include "auxiliary.h"
#include "kmc_base.h"

class MdMdkmcBindUnbind : public KMCBase {
  protected:
    double eps_eff_;
    double on_rate_;

    int nsimples_;
    std::vector<Simple*>* simples_;

    void Bind();
    void Unbind(SpeciesBase* spec);
    void FinishKMC();

  public:
    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking, int ikmc, YAML::Node &node, long seed);
    virtual void RunKMC(SpeciesBase *spec1, SpeciesBase *spec2);
    virtual void Print();
};

#endif
