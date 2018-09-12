#ifndef _SIMCORE_INTERACTION_H_
#define _SIMCORE_INTERACTION_H_
#include <tuple>

class Object;

typedef std::pair<int,int> oid_pair;
typedef std::pair<int,int> mesh_pair;

class Interaction {
  public:
    Interaction() {}
    oid_pair oids;
    mesh_pair mids;
    bool boundary = false; // true if boundary interaction
    double force[3] = {0}, // force acting on obj1 due to obj2
           t1[3] = {0}, // torque acting on obj1
           t2[3] = {0}, // torque acting on obj2
           dr[3] = {0}, // vector from obj1 to obj2
           contact1[3] = {0}, // vector from obj1 COM along obj1 to intersection with dr
           contact2[3] = {0}, // vector from obj2 COM along obj2 to intersection with dr
           buffer_mag = 0, // sum of object radii
           buffer_mag2 = 0, // " " " " squared
           dr_mag2 = -1, // magnitude of dr vector squared
           stress[9] = {0}, // stress tensor for calculating pressure
           pote = 0; // potential energy
};

typedef std::pair<Object*,Object*> interactor_pair;
typedef std::pair< interactor_pair , Interaction > pair_interaction;

#endif
