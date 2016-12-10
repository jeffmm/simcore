#ifndef _SIMCORE_INTERACTION_H_
#define _SIMCORE_INTERACTION_H_

class Interaction {
  public:
    Interaction() {}
    double force[3] = {0}, // force acting on obj1 due to obj2
           t1[3] = {0}, // torque acting on obj1
           t2[3] = {0}, // torque acting on obj2
           dr[3] = {0}, // vector from obj1 to obj2
           contact1[3] = {0}, // vector from obj1 COM along obj 1to intersection with dr
           contact2[3] = {0}, // vector from obj2 COM along obj2 to intersection with dr
           buffer_mag = 0, // sum of object radii
           buffer_mag2 = 0, // " " " " squared
           dr_mag2 = 0, // magnitude of dr vector squared
           virial = 0, // virial for calculating pressure
           pote = 0; // potential energy
};

#endif
