#ifndef _CYTOSCORE_SIMPLE_H_
#define _CYTOSCORE_SIMPLE_H_

#include "object.h"
#include "auxiliary.h"

class Simple : public Object {
  public:
    Simple(system_parameters *params, space_struct *space, long seed) :
      Object(params, space, seed) {
        is_simple_ = true;
    }
    virtual ~Simple() {}
    Simple(const Simple& that) : Object(that) {}
    Simple& operator=(Simple const& that) {
      Object::operator=(that);
      return *this;
    }
};

#endif // _CYTOSCORE_SIMPLE_H_
