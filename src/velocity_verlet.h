#ifndef _SIMCORE_VELOCITY_VERLET_H_
#define _SIMCORE_VELOCITY_VERLET_H_

#include "integrator.h"

// Velocity verlet integrator
// x(t+dt) = x(t) + v(t)*dt + 0.5*F(t)*dt^2/m
// calculate F(t+dt)
// v(t+dt) = v(t) + 0.5*(F(t) + F(t+dt))*dt/m

template <typename T>
class VelocityVerlet : public Integrator<T> {
  protected:
    using Integrator<T>::objects_;
    using Integrator<T>::params_;
//    using Integrator<T>::forces_;
  public:
    VelocityVerlet() : Integrator<T>() {}
    VelocityVerlet(system_parameters *params, std::vector<T> *objects, long seed) :
      Integrator<T>(params, objects, seed) {}
    void CalculateAllForces();
    void UpdateEnergies();
    void InitForces();
    void Integrate();
};

#endif // _SIMCORE_VELOCITY_VERLET_H_
