// Brute force calculation of forces
// O(N^2) at best

#ifndef _CYTOSCORE_FORCE_BRUTE_H_
#define _CYTOSCORE_FORCE_BRUTE_H_

#include "force_base.h"

class ForceBrute : public ForceBase {
  public:

    ForceBrute() {}
    virtual ~ForceBrute() {}

    // Init is the same
    // LoadSimples is the same

    virtual void InitMP();
    virtual void Finalize();
    virtual void UpdateScheme();
    virtual void Interact();
};

#endif
