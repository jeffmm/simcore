#ifndef _SIMCORE_BOND_H_
#define _SIMCORE_BOND_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

typedef enum {
  SHRINK = 0,
  GROW = 1,
  PAUSE
} poly_state_t;

class Bond : public Rigid {
  public:
    Bond(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Rigid(params, space, seed, sid) {}
    ~Bond() {}
    Bond(const Bond& that) : Rigid(that) {}
    Bond& operator=(Bond const& that) {
      Rigid::operator=(that); return *this;
    }
    void Init();

};

#endif // _SIMCORE_BOND_H_
