#ifndef _CYTOSCORE_BEAD_H_
#define _CYTOSCORE_BEAD_H_

#include "object.h"
#include "auxiliary.h"

class Bead : public Simple {
  public:
    Bead(system_parameters *params, space_struct * space, long seed, SID sid) : Simple(params, space, seed, sid) {}
    ~Bead() {}
    Bead(const Bead& that) : Simple(that) {}
    Bead& operator=(Bead const& that) {Simple::operator=(that); return *this;} 
    void InsertRandom() {
      generate_random_unit_vector(n_dim_, position_, rng_.r);
      double mag = gsl_rng_uniform_pos(rng_.r) * space_->radius;
      for (int i=0; i<n_dim_; ++i) {
        position_[i] *= mag;
        orientation_[i] = 0;
      }
    }
};

#endif // _CYTOSCORE_BEAD_H_
