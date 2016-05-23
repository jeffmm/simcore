#ifndef _CYTOSCORE_COMPOSITE_H_
#define _CYTOSCORE_COMPOSITE_H_

#include "auxiliary.h"
#include "object.h"
#include "parameters.h"

template <typename T>
class Composite : public Object {
  protected:
    std::vector<T> elements_;
  public:
    Composite(system_parameters *params, space_struct *space, long seed) : Object(params, space, seed) {
      is_simple_=false;
    } 
    //Destructor
    virtual ~Composite() {}
    //Copy constructor
    Composite(const Composite& that) : Object(that) {
      elements_=that.elements_;
    }
    //Assignment constructor
    Composite& operator=(Composite const& that) {
      Object::operator=(that);
      elements_=that.elements_;
      return *this;
    }
    virtual void ZeroForce() {
      memset(force_, 0, sizeof(force_));
      for (auto it=elements_.begin(); it!=elements_.end(); ++it)
        it->ZeroForce();
    }
};
#endif // _CYTOSCORE_COMPOSITE_H_
