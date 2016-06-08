// Microcell force computation

#ifndef _CYTOSCORE_FORCE_MICROCELL_H_
#define _CYTOSCORE_FORCE_MICROCELL_H_

#include "force_base.h"
#include "microcell_list.h"

class ForceMicrocell : public ForceBase {
  public:

    ForceMicrocell() {}
    virtual ~ForceMicrocell() {}

    virtual void Init(space_struct* pSpace, double pSkin) {
        // Override this to call base class, then initmp
        ForceBase::Init(pSpace, pSkin);
        InitMP();
    }

    virtual void InitMP();
    virtual void UpdateScheme();
    virtual void Interact();

  protected:

    MicrocellList microcell_list_;
};

#endif
