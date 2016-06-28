// Microcell force computation

#ifndef _SIMCORE_FORCE_MICROCELL_H_
#define _SIMCORE_FORCE_MICROCELL_H_

#include "force_base.h"
#include "microcell_list.h"

class ForceMicrocell : public ForceBase {
  public:

    ForceMicrocell() {
        name_ = "ForceMicrocell";
    }
    virtual ~ForceMicrocell() {}

    // Override these functions, need special stuff
    virtual void Init(space_struct* pSpace, std::vector<SpeciesBase*> *pSpecies, std::vector<Simple*> *pSimples, double pSkin);
    virtual void LoadSimples();

    virtual void printSpecifics();
    virtual void dump();

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    MicrocellList microcell_list_;
};

#endif
