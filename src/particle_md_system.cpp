#include "particle_md_system.h"


void ParticleMDSystem::Integrate() {
  integrator_.Integrate();
}

void ParticleMDSystem::InitElements() {
  int n_particles = params_->n_particles;
  double mass = params_->particle_mass;
  for (int i=0; i<n_particles; ++i)
    AddElement();
  double *vec = new double [n_dim_];
  double *net_v = new double [n_dim_];
  for (int i=0; i<n_dim_; ++i)
    net_v[i] = 0.0;
  for (auto elem=elements_.begin(); elem!=elements_.end(); ++elem) {
    elem->SetMass(mass);
    // Set random coordinate
    space_->RandomCoordinate(vec);
    elem->SetPosition(vec);
    // Set uniformly random velocity
    for (int i=0; i<n_dim_; ++i)
      net_v[i] += vec[i] = gsl_rng_uniform_pos(rng_.r) - 0.5;
    elem->SetVelocity(vec);
  }
  // Enforce zero net linear momentum
  for (int i=0; i<n_dim_; ++i)
    net_v[i]/=n_particles;
  for (auto elem=elements_.begin(); elem!=elements_.end(); ++elem) {
    double *vel = elem->GetVelocity();
    for (int i=0; i<n_dim_; ++i)
      vec[i] = vel[i] - net_v[i];
    elem->SetVelocity(vec);
  }
  delete[] vec;
  delete[] net_v;
  ZeroForces();
}

