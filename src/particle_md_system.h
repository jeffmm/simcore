#ifndef _CYTOSCORE_PARTICLE_MD_SYSTEM_H_
#define _CYTOSCORE_PARTICLE_MD_SYSTEM_H_

#include "object_system.h"
#include "particle_md.h"
#include "velocity_verlet.h"

class ParticleMDSystem : 
  public ObjectSystemBase<ParticleMD, VelocityVerlet> {

  public:
    ParticleMDSystem() : ObjectSystemBase< ParticleMD, VelocityVerlet>() {}
    ParticleMDSystem(system_parameters *params, SpaceProperties *space, long seed) : ObjectSystemBase< ParticleMD, VelocityVerlet>(params, space, seed) {
  SetTypeID(100);
  simple_ = true;
}
    ~ParticleMDSystem() {}
    void Integrate();
    void InitElements();
};

#endif // _CYTOSCORE_PARTICLE_MD_SYSTEM_H_
