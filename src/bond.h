#ifndef _SIMCORE_BOND_H_
#define _SIMCORE_BOND_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

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
    void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_BOND_H_
