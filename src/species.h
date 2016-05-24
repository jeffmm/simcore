#ifndef _CYTOSCORE_SPECIES_H_
#define _CYTOSCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"

class SpeciesBase {
  private:
    unsigned int sid_;
    static unsigned int next_sid_;
  protected:
    int n_members_;
    system_parameters *params_;
    space_struct *space_;
    rng_properties rng_;
  public:
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed) {
      sid_=++next_sid_;
      n_members_ = n_members;
      params_ = params;
      space_ = space;
      rng_.init(seed);
    }
    virtual ~SpeciesBase() {rng_.clear();}
    SpeciesBase(const SpeciesBase& that) {
      next_sid_=that.next_sid_;
      sid_=that.sid_;
      params_=that.params_;
      space_=that.space_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    SpeciesBase& operator=(SpeciesBase const& that) {
      next_sid_=that.next_sid_;
      sid_=that.sid_;
      params_=that.params_;
      space_=that.space_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    virtual void UpdatePositions() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init() {}
    virtual void ReInit(unsigned int const cid) {}
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> sim;
      return sim;
    }
    virtual double GetTotalEnergy() {return 0;}
    unsigned int const GetSID() {return sid_;}
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
    // Virtual functions
    virtual void Init() {
      for (int i=0; i<n_members_; ++i) {
        T * member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
        member->Init();
        members_.push_back(member);
      }
    }

    virtual void ReInit(unsigned int const cid) {
      for (auto it=members_.begin(); it!=members_.end();) {
        if ((*it)->GetCID() == cid) {
          delete (*it);
          members_.erase(it);
          T * member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
          member->Init();
          members_.push_back(member);
          break;
        }
        else {
          ++it;
        }
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
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> simples;
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        std::vector<Simple*> sim_vec = (*it)->GetSimples();
        simples.insert(simples.end(), sim_vec.begin(), sim_vec.end());
      }
      return simples;
    }
};

#endif // _CYTOSCORE_SPECIES_H_
