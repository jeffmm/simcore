#ifndef _SIMCORE_SITE_H_
#define _SIMCORE_SITE_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

class Site : public Simple {
  public:
    Site(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {}
    ~Site() {}
    Site(const Site& that) : Simple(that) {}
    Site& operator=(Site const& that) {
      Simple::operator=(that); return *this;
    }
    void Init();
    void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_SITE_H_
