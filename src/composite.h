#ifndef _CYTOSCORE_COMPOSITE_H_
#define _CYTOSCORE_COMPOSITE_H_

#include "auxiliary.h"
#include "simple.h"
#include "parameters.h"

class CompositeBase {
  protected:
    int n_dim_;
    int n_periodic_;
    double delta_;
    rng_properties rng_;
    graph_struct g_;
    space_struct *space_;
    double length_;
    double diameter_;
    double position_[3];
    double scaled_position_[3];
    double orientation_[3];
    double force_[3];

  public:
    CompositeBase(int n_dim, double delta, long seed) {
      n_dim_=n_dim;
      delta_ = delta;
      memset(position_,0,sizeof(position_));
      memset(scaled_position_,0,sizeof(scaled_position_));
      memset(orientation_,0,sizeof(orientation_));
      memset(force_,0,sizeof(force_));
      rng_.init(seed);
      // Set defaults
      diameter_ = 1;
      length_ = 0;
    }
    CompositeBase(const CompositeBase& that) {
      n_dim_ = that.n_dim_;
      delta_ = that.delta_;
      length_ = that.length_;
      diameter_ = that.diameter_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(scaled_position_,that.scaled_position_,sizeof(scaled_position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(force_,that.force_,sizeof(force_));
      rng_.init(gsl_rng_get(that.rng_.r));
      g_ = that.g_;
    }
    CompositeBase& operator=(CompositeBase const& that) {
      n_dim_ = that.n_dim_;
      delta_ = that.delta_;
      memcpy(position_,that.position_,sizeof(position_));
      memcpy(scaled_position_,that.scaled_position_,sizeof(scaled_position_));
      memcpy(orientation_,that.orientation_,sizeof(orientation_));
      memcpy(force_,that.force_,sizeof(force_));
      length_=that.length_;
      diameter_=that.diameter_;
      rng_.init(gsl_rng_get(that.rng_.r));
      g_ = that.g_;
      return *this;
    }
    double const GetLength() {return length_;}
    double const GetDiameter() {return diameter_;}
    double const * const GetPosition() {return position_;}
    double const * const GetScaledPosition() {return scaled_position_;}
    double const * const GetOrientation() {return orientation_;}
    void SetLength(double const length) {length_ = length;}
    void SetDiameter(double const diameter) {diameter_ = diameter;}
    void SetPosition(const double *const pos) {
      memcpy(position_,pos,n_dim_*sizeof(double));
    }
    void SetScaledPosition(const double *const scaled_pos) {
      memcpy(scaled_position_,scaled_pos,n_dim_*sizeof(double));
    }
    void SetOrientation(const double *const u) {
      memcpy(orientation_,u,n_dim_*sizeof(double));
    }
    virtual void UpdatePosition() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init(system_parameters *params, space_struct *space) {}
    virtual ~CompositeBase() {rng_.clear();}
};

template <typename T>
class Composite : public CompositeBase {
  protected:
    std::vector<T> elements_;
  public:
    Composite(int n_dim, double delta, long seed) : CompositeBase(n_dim, delta, seed) {} 
    //Destructor
    virtual ~Composite() {}
    //Copy constructor
    Composite(const Composite& that) : CompositeBase(that) {
      elements_=that.elements_;
    }
    //Assignment constructor
    Composite& operator=(Composite const& that) {
      CompositeBase::operator=(that);
      elements_=that.elements_;
      return *this;
    }
    virtual void ZeroForces() {
      for (auto it=elements_.begin(); it!=elements_.end(); ++it)
        it->ZeroForce();
    }
    virtual void UpdatePeriodic() {
      if (space_->n_periodic == 0)
        return;
      for (auto it=elements_.begin(); it!=elements_.end(); ++it) {
        double r[3],s[3];
        memcpy(r,it->GetPosition(),n_dim_*sizeof(double));
        memcpy(s,it->GetScaledPosition(),n_dim_*sizeof(double));
        periodic_boundary_conditions(space_->n_periodic, space_->unit_cell, space_->unit_cell_inv, r, s);
        it->SetPosition(r);
        it->SetScaledPosition(s);
      }
    }
};
#endif // _CYTOSCORE_COMPOSITE_H_
