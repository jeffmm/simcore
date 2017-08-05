#ifndef _SIMCORE_PARTICLE_H_
#define _SIMCORE_PARTICLE_H_

#include "object.h"

class Particle : public Object {
  public:
    Particle(int n_dim) : Object(n_dim);
};

Particle::Particle(int n_dim) : Object(n_dim) {
  SetTypeID(100);
  simple_ = true;
}

#endif // _SIMCORE_PARTICLE_H_
