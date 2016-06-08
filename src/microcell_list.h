// Microcell list!

#ifndef _CYTOSCORE_MICROCELL_LIST_H_
#define _CYTOSCORE_MICROCELL_LIST_H_

#include "force_substructure_base.h"

class MicrocellList : public ForceSubstructureBase {
  public:

    MicrocellList() {}
    virtual ~MicrocellList() {}

    // INIT SHOULD be the same
    // Load Simples should be the same
    virtual void CreateSubstructure(double pRcut);
};

#endif
