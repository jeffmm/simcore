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
  private:
    double graph_diameter_;
  public:
    Bond(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Rigid(params, space, seed, sid) {
      graph_diameter_ = params->graph_diameter;
    }
    void Init();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_BOND_H_
