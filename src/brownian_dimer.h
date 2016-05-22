#ifndef _CYTOSCORE_BROWNIAN_DIMER_H_
#define _CYTOSCORE_BROWNIAN_DIMER_H_

#include "auxiliary.h"
#include "parameters.h"
#include <yaml-cpp/yaml.h>
#include "composite.h"
#include "bead.h"

class BrownianDimer : public Composite<Bead> {
  private:
    /* Unique to BrownianDimer */
    double eq_length_;
    double k_spring_;

    /* Functions unique to Brownian Dimer */
    void Integrate();
    void InternalForces();
    void KickBeads();
    void UpdateOrientation();
    void InitRandom(double sys_radius); // Temporary, get space to do this
    graph_struct g2_; // Drawing two spheros, one for each bead;

  public:
    //Constructor
    BrownianDimer(int n_dim, double delta, long seed) : Composite(n_dim, delta, seed) {
      for (int i=0; i<2; ++i) {
        Bead b(n_dim_, delta_, gsl_rng_get(rng_.r));
        elements_.push_back(b);
      }
      // Set defaults
      diameter_ = 1;
      length_ = 5;
      eq_length_ = 5;
      k_spring_ = 10;
    }
    //Destructor
    ~BrownianDimer() {}
    //Copy constructor
    BrownianDimer(const BrownianDimer& that) : Composite(that) {
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
      space_=that.space_;
    }
    //Assignment constructor
    BrownianDimer& operator=(BrownianDimer const& that) {
      Composite::operator=(that);
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
      space_=that.space_;
      return *this;
    };
    /* Define virtual functions */
    void Init(system_parameters *params, space_struct *space) {
      diameter_ = params->dimer_diameter;
      length_ = params->dimer_length;
      eq_length_ = params->dimer_eq_length;
      k_spring_ = params->dimer_k_spring;
      space_=space;
      InitRandom(space_->radius);
    }
    void UpdatePosition();
    void Draw(std::vector<graph_struct*> * graph_array);

    /* Functions unique to Brownian Dimer */
    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
};

#endif // _CYTOSCORE_BROWNIAN_DIMER_H_
