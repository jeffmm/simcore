#include "integrator_manager.h"

IntegratorManager::IntegratorManager() {}
IntegratorManager::~IntegratorManager() {}

void IntegratorManager::Init(system_parameters *params, std::vector<Objects*> *objects, long seed) {
  GetMethod(params, objects, seed);
}

void IntegratorManager::GetMethod(system_parameters *params, std::vector<Objects*> *objects, long seed) {
  if (params_->integrator.compare("velocity_verlet")) {
    method = new velocity_verlet(params, objects, seed);
  }
  else {
    error_exit("ERROR: Integrator not recognized.\n");
  }
}

void IntegratorManager::Clear() {
  delete method_;
}

void IntegratorManager::Integrate(std::vector<Object*> *objects) {
  method->Integrate(objects);
}

