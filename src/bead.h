#ifndef _CYTOSCORE_BEAD_H_
#define _CYTOSCORE_BEAD_H_

#include "auxiliary.h"

class Bead {
  private:
    int n_dim_;
    double delta_;
    double mass_;
    double position_[3];
    double orientation_[3];
    double prev_position_[3];
    double velocity_[3];
    double force_[3];
    double diameter_;
    double length_;
    graph_struct g_;
    rng_properties rng_;
  public:
    Bead(int n_dim, double delta, long seed) {
      n_dim_=n_dim; 
      delta_ = delta;
      memset(position_,0,sizeof(position_));
      memset(orientation_,0,sizeof(orientation_));
      memset(prev_position_,0,sizeof(prev_position_));
      memset(velocity_,0,sizeof(velocity_));
      memset(force_,0,sizeof(force_));
      mass_ = 1;
      diameter_ = 1;
      length_ = 0;
      rng_.init(seed);
    }
    ~Bead() {rng_.clear();}
    Bead(const Bead& that) {
      n_dim_=that.n_dim_;
      delta_ = that.delta_;
      mass_ = that.mass_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(prev_position_,that.prev_position_,sizeof(prev_position_));
      memcpy(velocity_,that.velocity_,sizeof(velocity_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_ = that.diameter_;
      length_ = that.length_;
      g_ = that.g_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    Bead& operator=(Bead const& that) {
      n_dim_=that.n_dim_;
      delta_ = that.delta_;
      mass_ = that.mass_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(prev_position_,that.prev_position_,sizeof(prev_position_));
      memcpy(velocity_,that.velocity_,sizeof(velocity_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_ = that.diameter_;
      length_ = that.length_;
      g_ = that.g_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }

    void SetNDim(int n_dim) {n_dim_ = n_dim;}
    void SetPosition(const double *const pos) {
      memcpy(position_,pos,n_dim_*sizeof(double));
    }
    void SetOrientation(const double *const u) {
      memcpy(orientation_,u,n_dim_*sizeof(double));
    }
    void SetPrevPosition() {
      memcpy(prev_position_,position_,sizeof(position_));
    }
    void SetMass(const double m) {mass_ = m;}
    void SetDiameter(const double r) {diameter_ = r;}    
    void SetLength(const double l) {length_ = l;}
    void ZeroForce() {memset(force_,0,sizeof(force_));}
    void AddForce(const double *const f) {
      for (int i=0; i<n_dim_; ++i)
        force_[i]+=f[i];
    }
    double const * const GetForce() {return force_;}
    double const * const GetPosition() {return position_;}
    double const * const GetPrevPosition() {return prev_position_;}
    double const * const GetOrientation() {return orientation_;}
    double const * const GetVelocity() {return velocity_;}
    double const GetMass() {return mass_;}
    double const GetLength() {return length_;}
    double const GetDiameter() {return diameter_;}
    void Draw(std::vector<graph_struct*> * graph_array) {
      memcpy(g_.r,position_,sizeof(position_));
      memcpy(g_.u,orientation_,sizeof(orientation_));
      g_.length = length_;
      g_.diameter = diameter_;
      graph_array->push_back(&g_);
    }
};

#endif // _CYTOSCORE_BEAD_H_
