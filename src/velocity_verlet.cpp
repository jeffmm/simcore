
#include "velocity_verlet.h"

template <typename T>
void VelocityVerlet<T>::InitForces() {
  CalculateAllForces(); // Get F(t)
}

template <typename T>
void VelocityVerlet<T>::Integrate() {
  double *pos, *vel, *force;
  for (auto ob=objects_->begin(); ob!=objects_->end(); ++ob) {
    pos = ob->GetPosition();
    vel = ob->GetVelocity();
    force = ob->GetForce();
    for (int i=0; i<params_->n_dim; ++i) {
      // Update position using current force
      pos[i] = pos[i] + params_->delta * vel[i] + 0.5 * force[i] * SQR(params_->delta) / ob->mass;
      // Update velocity halfstep using current force
      vel[i] = vel[i] + 0.5 * force[i] * params_->delta / ob->mass;
    }
  }
  // Calculate new forces and finish velocity update
  CalculateAllForces();
  for (auto ob=objects_->begin(); ob!=objects_->end(); ++ob) {
    vel = ob->GetVelocity();
    for (int i=0; i<params_->n_dim; ++i)
      vel[i] = vel[i] + 0.5 * force[i] * params_->delta / ob->mass;
  }
  // Update particle energies
  UpdateEnergies();
}

template <typename T>
void VelocityVerlet<T>::CalculateAllForces() {
  //forces_.GetInteractions(objects_);
  for (auto ob=objects_->begin(); ob!=objects_->end(); ++ob) {
    ob->ApplyInteractions();
  }
}

template <typename T>
void VelocityVerlet<T>::UpdateEnergies() {
  for (auto ob=objects_->begin(); ob!=objects_->end(); ++ob) {
    ob->CalculateEnergy();
  }
}

