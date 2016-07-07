#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"
#include "potential_base.h"

class SpeciesBase {
  private:
    SID sid_;
  protected:
    int n_members_;
    bool is_kmc_;
    system_parameters *params_;
    space_struct *space_;
    rng_properties rng_;
    void SetSID(SID sid) {sid_=sid;}
    std::vector<potential_pair> potentials_;
    //virtual void InitPotentials(system_parameters *params) {}
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
      is_kmc_ = false;
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
      is_kmc_ = that.is_kmc_;
      rng_.init(gsl_rng_get(that.rng_.r));
    }
    SpeciesBase& operator=(SpeciesBase const& that) {
      sid_=that.sid_;
      params_=that.params_;
      space_=that.space_;
      potentials_=that.potentials_;
      is_kmc_=that.is_kmc_;
      rng_.init(gsl_rng_get(that.rng_.r));
        return *this;
    }
    virtual void UpdatePositions() {}
    virtual void UpdatePositionsMP() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init() {}
    virtual void ReInit(unsigned int const cid) {}
    virtual double GetDrMax() {return 0.0;}
    virtual void ZeroDr() {};
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> sim;
      return sim;
    }
    virtual double GetKineticEnergy() {return 0;}
    virtual double GetPotentialEnergy() {return 0;}
    virtual double GetTotalEnergy() {return 0;}
    SID const GetSID() {return sid_;}
    bool IsKMC() {return is_kmc_;}
    virtual void PrepKMC() {}
    virtual void StepKMC() {}
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
    virtual void UpdatePositionsMP() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->UpdatePositionMP();
      }
    }
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> simples;
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        std::vector<Simple*> sim_vec = (*it)->GetSimples();
        simples.insert(simples.end(), sim_vec.begin(), sim_vec.end());
      }
      return simples;
    }
    virtual double GetKineticEnergy() {
      double ke=0;
      for (auto it=members_.begin(); it!=members_.end(); ++it)
        ke+=(*it)->GetKineticEnergy();
      return ke;
    }
    virtual double GetPotentialEnergy() {
      double pe=0;
      for (auto it=members_.begin(); it!=members_.end(); ++it)
        pe+=(*it)->GetPotentialEnergy();
      // The total potential energy is going to be half of the
      // potential energy felt by each particle. Potential energy is shared,
      // so I need to avoid double counting.
      return 0.5*pe;
    }
    virtual double GetTotalEnergy() {
      double ke = GetKineticEnergy();
      double pe = GetPotentialEnergy();
      return ke+pe;
    }
    virtual double GetDrMax() {
      double dr_mag2_max = 0.0;
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        double dr_mag2 = 0.0;
        double const * const dr = (*it)->GetDrTot();
        for (int i=0; i<space_->n_dim; ++i)
          dr_mag2 += dr[i]*dr[i];
        if (dr_mag2 > dr_mag2_max)
          dr_mag2_max = dr_mag2;
      }
      return dr_mag2_max;
    }
    virtual void ZeroDr() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->ZeroDrTot();
      }
    }
};

#endif // _SIMCORE_SPECIES_H_
