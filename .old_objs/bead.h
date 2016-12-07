#ifndef _SIMCORE_BEAD_H_
#define _SIMCORE_BEAD_H_

#include "object.h"
#include "auxiliary.h"

class Bead : public Simple {
  public:
    Bead(system_parameters *params, space_struct * space, long seed, SID sid) : Simple(params, space, seed, sid) {}
    ~Bead() {}
    Bead(const Bead& that) : Simple(that) {}
    Bead& operator=(Bead const& that) {Simple::operator=(that); return *this;} 
};

#endif // _SIMCORE_BEAD_H_
