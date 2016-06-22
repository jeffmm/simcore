#ifndef _SIMCORE_BOND_H_
#define _SIMCORE_BOND_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

class Bond : public Simple {
  public:
    Bond(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {}
    ~Bond() {}
    Bond(const Bond& that) : Simple(that) {}
    Bond& operator=(Bond const& that) {
      Simple::operator=(that); return *this;
    }
    void Init();
    void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_BOND_H_
