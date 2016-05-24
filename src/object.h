#ifndef _CYTOSCORE_OBJECT_H_
#define _CYTOSCORE_OBJECT_H_

#include "auxiliary.h"

class Object {

  private:
    unsigned int oid_;
    static unsigned int next_oid_;
    static unsigned int next_cid_;
        
  protected:
    unsigned int cid_;
    unsigned int sid_;
    int n_dim_;
    double position_[3],
           scaled_position_[3],
           prev_position_[3],
           orientation_[3],
           force_[3],
           delta_,
           diameter_,
           length_;
    bool is_simple_;
    space_struct *space_;
    graph_struct g_;
    rng_properties rng_;
    virtual void InsertRandom(double buffer);
    unsigned int const NextCID() {return ++next_cid_;}
  public:
    Object(system_parameters *params, space_struct *space, long seed, unsigned int const sid);
    Object(const Object& that);
    Object& operator=(Object const& that);

    virtual ~Object() {rng_.clear();}
    bool IsSimple() {return is_simple_;}
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
    void SetDiameter(double new_diameter) {diameter_ = new_diameter;}
    void SetLength(double new_length) {length_ = new_length;}
    void SetSpace(space_struct * space) {space_ = space;}
    void ZeroForce() {memset(force_,0,sizeof(force_));}
    void AddForce(const double *const f) {
      for (int i=0; i<n_dim_; ++i)
        force_[i]+=f[i];
    }
    double const * const GetPosition() {return position_;}
    double const * const GetScaledPosition() {return scaled_position_;}
    double const * const GetOrientation() {return orientation_;}
    double const * const GetForce() {return force_;}
    double const GetDiameter() {return diameter_;}
    double const GetLength() {return length_;}
    virtual void Init() {InsertRandom(length_+diameter_);}
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePeriodic();
    virtual void UpdatePosition() {};
    unsigned int const GetOID() {return oid_;}
    unsigned int const GetCID() {return cid_;}
    unsigned int const GetSID() {return sid_;}
};

class Simple : public Object {
  public:
    Simple(system_parameters *params, space_struct *space, long seed, unsigned int const sid) :
      Object(params, space, seed, sid) {
        SetCID(GetOID());
        is_simple_ = true;
    }
    virtual ~Simple() {}
    Simple(const Simple& that) : Object(that) {}
    Simple& operator=(Simple const& that) {
      Object::operator=(that);
      return *this;
    }
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> sim_vec;
      sim_vec.push_back(this);
      return sim_vec;
    }
    void SetCID(unsigned int const cid) {cid_=cid;}
};

template <typename T>
class Composite : public Object {
  protected:
    std::vector<T> elements_;
  public:
    Composite(system_parameters *params, space_struct *space, long seed, unsigned int sid) : Object(params, space, seed, sid) {
      cid_ = NextCID();
      is_simple_=false;
    } 
    //Destructor
    virtual ~Composite() {}
    //Copy constructor
    Composite(const Composite& that) : Object(that) {
      elements_=that.elements_;
    }
    //Assignment constructor
    Composite& operator=(Composite const& that) {
      Object::operator=(that);
      elements_=that.elements_;
      return *this;
    }
    virtual void ZeroForce() {
      memset(force_, 0, sizeof(force_));
      for (auto it=elements_.begin(); it!=elements_.end(); ++it)
        it->ZeroForce();
    }
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> sim_vec;
      for (auto it=elements_.begin(); it!=elements_.end(); ++it)
        sim_vec.push_back(&(*it));
      return sim_vec;
    }
};


#endif // _CYTOSCORE_OBJECT_H_
