#ifndef _SIMCORE_INTEGRATOR_MANAGER_H_
#define _SIMCORE_INTEGRATOR_MANAGER_H_

#include "integrator.h"

class IntegratorManager {

  private:
    Integrator *method_;
    rng_properties rng_;
    system_parameters *params_;
  public:
    IntegratorManager();
    ~IntegratorManager();
    void Init(system_parameters *params, std::vector<Objects*> *objects, long seed);
    void GetMethod(system_parameters *params, std::vector<Objects*> *objects, long seed);
    void Clear();
};


#endif // _SIMCORE_INTEGRATOR_MANAGER_H_
