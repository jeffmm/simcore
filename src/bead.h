#ifndef _CYTOSCORE_BEAD_H_
#define _CYTOSCORE_BEAD_H_

#include "simple.h"
#include "auxiliary.h"

class Bead : public Simple {
  public:
    Bead(int n_dim, double delta, long seed) : Simple(n_dim, delta, seed) {}
    ~Bead() {}
    Bead(const Bead& that) : Simple(that) {}
    Bead& operator=(Bead const& that) {Simple::operator=(that); return *this;} 
    void Draw(std::vector<graph_struct*> * graph_array) {
      memcpy(g_.r,position_,sizeof(position_));
      memcpy(g_.u,orientation_,sizeof(orientation_));
      g_.length = length_;
      g_.diameter = diameter_;
      graph_array->push_back(&g_);
    }
    void KickBead() {
      for (int i=0; i<n_dim_; ++i) {
        double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
        force_[i] += kick*diffusion_;
      }
    }
};

#endif // _CYTOSCORE_BEAD_H_
