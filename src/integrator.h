#ifndef _SIMCORE_INTEGRATOR_H_
#define _SIMCORE_INTEGRATOR_H_

#include "auxiliary.h"
#include "object.h"

template <class T>
class Integrator {
  protected:
    system_parameters *params_;
    rng_properties rng_;
    std::vector<T> *objects_;
  public:
    Integrator() {}
    virtual ~Integrator() { rng_.clear(); }
    Integrator(system_parameters *params, std::vector<T> *objects, long seed) { Init(params, objects, seed); }    
    virtual void Init(system_parameters *params, std::vector<T> *objects, long seed) {
      params_ = params;
      objects_ = objects;
      rng_.init(seed);
      //InitForces();
    }
    virtual void InitForces() {}
    virtual void Integrate() {}
};



#endif // _SIMCORE_INTEGRATOR_H_
