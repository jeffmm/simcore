#ifndef _CYTOSCORE_OBJECT_SYSTEM_H_
#define _CYTOSCORE_OBJECT_SYSTEM_H_

#include "auxiliary.h"
#include "space.h"
#include "object.h"
class ObjectSystem {

  protected:
    int n_dim_;
    bool simple_;
    unsigned int tid_;
    int n_elements_;
    system_parameters *params_;
    SpaceProperties *space_;
    rng_properties rng_;
    virtual void SetTypeID(unsigned int new_tid);

  public:
    ObjectSystem() {}
    virtual ~ObjectSystem() {}
    bool IsSimple();
    const unsigned int GetTypeID() const;
    virtual void ZeroForces() {}
    virtual void Init() {};
    virtual void Integrate() {}
    virtual void AddElement() {}
    int GetNElements() {return n_elements_;}
    virtual void GetGraphics(graph_struct *g_struct) {}
    virtual void InitElements() {}
};

template <typename T, template<typename> class V>
class ObjectSystemBase : public ObjectSystem {
  protected:
    std::vector<T> elements_;
    V<T> integrator_;
  public:
    ObjectSystemBase() {}
    ObjectSystemBase(system_parameters *params, SpaceProperties *space, long seed) { 
      Init(params, space, seed);
    }
    virtual ~ObjectSystemBase();
    void Init(system_parameters *params, SpaceProperties *space, long seed);
    void ZeroForces();
    void AddElement();
    void GetGraphics(graph_struct *g_struct);
};


#endif // _CYTOSCORE_OBJECT_SYSTEM_H_
