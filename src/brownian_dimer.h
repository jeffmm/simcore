#ifndef _CYTOSCORE_BROWNIAN_DIMER_H_
#define _CYTOSCORE_BROWNIAN_DIMER_H_

#include "auxiliary.h"
#include "parameters.h"
#include <yaml-cpp/yaml.h>
#include "bead.h"

class BrownianDimer {
  private:
    /* Used in every composite object */
    int n_dim_;
    double delta_;
    rng_properties rng_;
    std::vector<Bead> beads_; // Generalize to simples
    graph_struct g_;
    double length_;
    double diameter_;
    double com_[3];
    double orientation_[3];

    /* Unique to BrownianDimer */
    double eq_length_;
    double k_spring_;
    double force_[3];
    double diffusion_;

    /* Functions used in every composite object */
    void ZeroForces();
    void Integrate();
    void InternalForces();

    /* Functions unique to Brownian Dimer */
    void KickBeads();
    void UpdateOrientation();
    void InitRandom(double sys_radius);

  public:
    //Constructor
    BrownianDimer(int n_dim, double delta, long seed) {
      n_dim_=n_dim;
      delta_ = delta;
      memset(com_,0,sizeof(com_));
      memset(orientation_,0,sizeof(orientation_));
      memset(force_,0,sizeof(force_));
      rng_.init(seed);
      for (int i=0; i<2; ++i) {
        Bead bead(n_dim_, delta_, gsl_rng_get(rng_.r));
        beads_.push_back(bead);
      }
      // Set defaults
      diameter_ = 1;
      length_ = 5;
      eq_length_ = 5;
      k_spring_ = 10;
    }
    //Destructor
    ~BrownianDimer() {rng_.clear();}
    //Copy constructor
    BrownianDimer(const BrownianDimer& that) {
      n_dim_ = that.n_dim_;
      delta_ = that.delta_;
      eq_length_ = that.eq_length_;
      k_spring_ = that.k_spring_;
      memcpy(com_,that.com_,sizeof(com_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_=that.diameter_;
      diffusion_=that.diffusion_;
      rng_.init(gsl_rng_get(that.rng_.r));
      beads_=that.beads_;
      g_ = that.g_;
    };
    //Assignment constructor
    BrownianDimer& operator=(BrownianDimer const& that) {
      n_dim_ = that.n_dim_;
      delta_ = that.delta_;
      eq_length_ = that.eq_length_;
      k_spring_ = that.k_spring_;
      memcpy(com_,that.com_,sizeof(com_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_=that.diameter_;
      diffusion_=that.diffusion_;
      rng_.init(gsl_rng_get(that.rng_.r));
      beads_=that.beads_;
      g_ = that.g_;
      return *this;
    };
    /* Initialization of parameters */
    void Init(system_parameters *params) {
      diameter_ = params->dimer_diameter;
      length_ = params->dimer_length;
      eq_length_ = params->dimer_eq_length;
      k_spring_ = params->dimer_k_spring;
      InitRandom(params->system_radius);
    }

    /* Used in every composite object */
    void UpdatePosition();
    void Draw(std::vector<graph_struct*> * graph_array);
    double const GetLength() {return length_;}
    double const GetDiameter() {return diameter_;}
    double const * const GetCOM() {return com_;}
    double const * const GetOrientation() {return orientation_;}
    void SetLength(double const length) {length_ = length;}
    void SetDiameter(double const diameter) {diameter_ = diameter;}

    /* Unique to Brownian Dimer */
    void InitRest();
    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
};

#endif // _CYTOSCORE_BROWNIAN_DIMER_H_
