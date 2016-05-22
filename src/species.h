#ifndef _CYTOSCORE_SPECIES_H_
#define _CYTOSCORE_SPECIES_H_

#include "auxiliary.h"
#include "parameters.h"

class SpeciesBase {
  protected:
    int n_members_;
    system_parameters *params_;
    space_struct *space_;
    rng_properties rng_;
  public:
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed) {
      n_members_ = n_members;
      params_ = params;
      space_ = space;
      rng_.init(seed);
    }
    virtual ~SpeciesBase() {rng_.clear();}
    SpeciesBase(const SpeciesBase& that) {
      params_=that.params_;
      space_=that.space_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    SpeciesBase& operator=(SpeciesBase const& that) {
      params_=that.params_;
      space_=that.space_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    virtual void UpdatePositions() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init() {}
};

template <typename T>
class Species : public SpeciesBase {
  protected:
    std::vector<T*> members_;
  public:
    Species(int n_members, system_parameters *params, space_struct *space, long seed) : SpeciesBase(n_members, params, space, seed) {}
    //Destructor
    virtual ~Species() {
      for (auto it=members_.begin(); it!=members_.end(); ++it)
        delete (*it);
    }
    //Copy constructor
    Species(const Species& that) : SpeciesBase(that) {
      members_=that.members_;
    }
    //Assignment constructor
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      members_=that.members_;
      return *this;
    }
    virtual void Init() {
      for (int i=0; i<n_members_; ++i) {
        T * member = new T(space_->n_dim, params_->delta, gsl_rng_get(rng_.r));
        member->Init(params_, space_);
        members_.push_back(member);
      }
    }
    virtual void Draw(std::vector<graph_struct*> * graph_array) {
      for (auto it=members_.begin(); it!=members_.end(); ++it)
        (*it)->Draw(graph_array);
    }
    virtual void UpdatePositions() {
      for (auto it=members_.begin(); it!=members_.end(); ++it)
        (*it)->UpdatePosition();
    }
};
#endif // _CYTOSCORE_SPECIES_H_
