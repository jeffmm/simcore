#ifndef _CYTOSCORE_SIMPLE_H_
#define _CYTOSCORE_SIMPLE_H_

#include "auxiliary.h"

class Simple {

  private:
    unsigned int oid_;
    static unsigned int next_oid_;
        
  protected:
    int n_dim_;
    double position_[3],
           scaled_position_[3],
           prev_position_[3],
           orientation_[3],
           velocity_[3],
           force_[3],
           delta_,
           diameter_,
           length_,
           mass_,
           diffusion_;
    graph_struct g_;
    rng_properties rng_;
    std::vector<interaction> interactions;

  public:
    Simple(int n_dim, double delta, long seed) {
      oid_ = ++next_oid_;
      n_dim_ = n_dim; 
      delta_ = delta;
      memset(position_,0,sizeof(position_));
      memset(scaled_position_,0,sizeof(scaled_position_));
      memset(orientation_,0,sizeof(orientation_));
      memset(prev_position_,0,sizeof(prev_position_));
      memset(velocity_,0,sizeof(velocity_));
      memset(force_,0,sizeof(force_));
      mass_ = 1;
      diameter_ = 1;
      length_ = 0;
      SetDiffusion();
      rng_.init(seed);
    }
    Simple(const Simple& that) {
      n_dim_=that.n_dim_;
      oid_ = that.oid_;
      next_oid_=that.next_oid_;
      delta_ = that.delta_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(scaled_position_,that.scaled_position_,sizeof(scaled_position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(prev_position_,that.prev_position_,sizeof(prev_position_));
      memcpy(velocity_,that.velocity_,sizeof(velocity_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_ = that.diameter_;
      length_ = that.length_;
      mass_ = that.mass_;
      diffusion_ = that.diffusion_;
      g_ = that.g_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    Simple& operator=(Simple const& that) {
      n_dim_=that.n_dim_;
      oid_=that.oid_;
      next_oid_=that.next_oid_;
      delta_ = that.delta_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(scaled_position_,that.scaled_position_,sizeof(scaled_position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(prev_position_,that.prev_position_,sizeof(prev_position_));
      memcpy(velocity_,that.velocity_,sizeof(velocity_));
      memcpy(force_,that.force_,sizeof(force_));
      diameter_ = that.diameter_;
      length_ = that.length_;
      mass_ = that.mass_;
      diffusion_ = that.diffusion_;
      g_ = that.g_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }

    virtual ~Simple() {rng_.clear();}
    const unsigned int GetOID() {return oid_;}
    void SetPosition(const double *const pos) {
      memcpy(position_,pos,n_dim_*sizeof(double));
    }
    void SetScaledPosition(const double *const scaled_pos) {
      memcpy(scaled_position_,scaled_pos,n_dim_*sizeof(double));
    }
    void SetOrientation(const double *const u) {
      memcpy(orientation_,u,n_dim_*sizeof(double));
    }
    void SetPrevPosition() {
      memcpy(prev_position_,position_,sizeof(position_));
    }
    void SetVelocity(const double *const vel) {
      memcpy(velocity_,vel,n_dim_*sizeof(double));
    }
    void SetDiameter(double new_diameter) {diameter_ = new_diameter;}
    void SetLength(double new_length) {length_ = new_length;}
    void SetMass(double new_mass) {mass_ = new_mass;}
    virtual void ApplyInteractions() {}
    void ZeroForce() {memset(force_,0,sizeof(force_));}
    void AddForce(const double *const f) {
      for (int i=0; i<n_dim_; ++i)
        force_[i]+=f[i];
    }
    void AddInteraction(interaction new_interaction) {
      interactions.push_back(new_interaction);
    }
    double const * const GetPosition() {return position_;}
    double const * const GetScaledPosition() {return scaled_position_;}
    double const * const GetOrientation() {return orientation_;}
    double const * const GetVelocity() {return velocity_;}
    double const * const GetForce() {return force_;}
    double const GetDiameter() {return diameter_;}
    double const GetLength() {return length_;}
    double const GetMass() {return mass_;}
    virtual void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
};

#endif // _CYTOSCORE_SIMPLE_H_
