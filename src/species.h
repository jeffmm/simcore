#ifndef _CYTOSCORE_SPECIES_H_
#define _CYTOSCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"
#include "potential_base.h"

class SpeciesBase {
  private:
    SID sid_;
  protected:
    int n_members_;
    system_parameters *params_;
    space_struct *space_;
    rng_properties rng_;
    void SetSID(SID sid) {sid_=sid;}
    std::vector<potential_pair> potentials_;
    virtual void InitPotentials() {}
    void AddPotential(SID sid1, SID sid2, PotentialBase * potential) {
      sid_pair sids = std::make_pair(sid1, sid2);
      potential_pair pot_pair = std::make_pair(sids,potential);
      potentials_.push_back(pot_pair);
    }
  public:
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed) {
      sid_ = SID::none;
      n_members_ = n_members;
      params_ = params;
      space_ = space;
      rng_.init(seed);
    }
    virtual ~SpeciesBase() {
      rng_.clear();
      for (auto it=potentials_.begin(); it!=potentials_.end(); ++it)
        delete it->second;
      potentials_.clear();
    }
    SpeciesBase(const SpeciesBase& that) {
      sid_=that.sid_;
      params_=that.params_;
      space_=that.space_;
      potentials_=that.potentials_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    SpeciesBase& operator=(SpeciesBase const& that) {
      sid_=that.sid_;
      params_=that.params_;
      space_=that.space_;
      potentials_=that.potentials_;
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
    SID const GetSID() {return sid_;}
    std::vector<potential_pair> GetPotentials() {return potentials_;}
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
