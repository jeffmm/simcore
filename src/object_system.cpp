
#include "object_system.h"


void ObjectSystemBase::SetTypeID(unsigned int new_tid) {
  tid_ = new_tid;
}

const unsigned int ObjectSystemBase::GetTypeID() const {
  return tid_;
}

bool ObjectSystemBase::IsSimple() {
  return simple_;
}

template <typename T, template<typename> class V>
void ObjectSystem<T,V>::Init(system_parameters *params, SpaceProperties *space, long seed) {
  n_elements_ = 0;
  n_dim_ = params->n_dim;
  params_ = params;
  space_ = space;
  rng_.init(seed);
  integrator_.Init(params, &elements_, gsl_rng_get(rng_.r));
}

template <typename T, template<typename> class V>
void ObjectSystem<T,V>::ZeroForces() {
  for (auto it=elements_.begin(); it!=elements_.end(); ++it)
    it->ZeroForce();
}

template <typename T, template<typename> class V>
void ObjectSystem<T,V>::AddElement() {
  n_elements_++;
  T* elem = new T(n_dim_);
  elements_.push_back(elem);
  delete elem;
}

template <typename T, template<typename> class V>
void ObjectSystem<T,V>::GetGraphics(graph_struct *g_struct) {
  double *r, *u;
  for (auto elem=elements_.begin(); elem!=elements_.end(); ++elem) {
    g_struct->l[g_struct->i_sphero] = elem->GetLength();
    g_struct->diam[g_struct->i_sphero] = 2.0*(elem->GetRadius());
    r = elem->GetPosition();
    u = elem->GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      g_struct->r[g_struct->i_sphero][i] = r[i];
      g_struct->u[g_struct->i_sphero][i] = u[i];
    }
    g_struct->i_sphero++;
  }
}

