#ifndef _SIMCORE_PARTICLE_MD_H_
#define _SIMCORE_PARTICLE_MD_H_

#include "object.h"

class ParticleMD : public Object {
  private:
    double energy_,
           mass_;
  public:
    ParticleMD(int n_dim) : Object(n_dim) {}
    void SetMass(double new_mass);
    double GetMass() const;
    void CalculateEnergy();
    double GetEnergy() const;
};

void ParticleMD::SetMass(double new_mass) {
  mass_ = new_mass;
}

void ParticleMD::CalculateEnergy() {
  double v_square;
  for (int i=0; i<n_dim_; ++i)
    v_square += SQR(velocity_[i]);
  energy_ = 0.5 * mass_ * v_square;
}

double ParticleMD::GetMass() const {
  return mass_;
}

double ParticleMD::GetEnergy() const {
  return energy_;
}

#endif // _SIMCORE_PARTICLE_MD_H_
