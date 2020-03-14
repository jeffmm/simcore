#ifndef _SIMCORE_INTERACTION_H_
#define _SIMCORE_INTERACTION_H_
#include "definitions.hpp"
#include <tuple>

class Object;

class Interaction {
public:
  Interaction() {}
  Interaction(Object *o1, Object *o2) : obj1(o1), obj2(o2) {}
  Interaction(Object *o1) : obj1(o1), boundary(true) {}
  Object *obj1 = nullptr;
  Object *obj2 = nullptr;
  bool boundary = false; // true if boundary interaction
  double force[3] = {0}; // force acting on obj1 due to obj2
  double t1[3] = {0};    // torque acting on obj1
  double t2[3] = {0};    // torque acting on obj2
  double dr[3] = {0};    // vector from obj1 to obj2
  double midpoint[3] = {0};
  // vector from obj1 COM along obj1 to intersection with dr
  double contact1[3] = {0};
  // vector from obj2 COM along obj2 to intersection with dr
  double contact2[3] = {0};
  double buffer_mag = 0;     // sum of object radii
  double buffer_mag2 = 0;    // " " " " squared
  double dr_mag2 = -1;       // magnitude of dr vector squared
  double stress[9] = {0};    // stress tensor for calculating pressure
  double pote = 0;           // potential energy
  double polar_order = 0;    // local polar order contribution
  double contact_number = 0; // contact number contribution
  bool no_interaction = false; // true for objects that do not interact
  bool pause_interaction = false; // no interaction this update
};

//typedef std::pair<int, Interaction*> object_interaction;
typedef std::pair<Interaction*, bool> object_interaction;

#endif
