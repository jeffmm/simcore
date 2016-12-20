#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"
#include <iterator>

class SpeciesBase {
  private:
    SID sid_;
  protected:
    bool posit_flag_;
    int n_members_,
        n_posit_;
    double delta_;
    system_parameters *params_;
    std::fstream oposit_file_;
    std::fstream iposit_file_;
    space_struct *space_;
    rng_properties rng_;
    void SetSID(SID sid) {sid_=sid;}

  public:
    SpeciesBase() {sid_ = SID::none;}
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed);
    virtual void UpdatePositions() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void InitConfig(system_parameters *params, space_struct *space, long seed);
    virtual void Init() {}
    virtual double GetDrMax() {return 0.0;}
    virtual void ZeroDr() {}
    virtual void ZeroForces() {}
    virtual std::vector<Simple*> GetSimples();
    virtual double GetKineticEnergy() {return 0;}
    virtual double GetPotentialEnergy() {return 0;}
    virtual double GetTotalEnergy() {return 0;}
    virtual void ScalePositions() {}
    SID const GetSID() {return sid_;}
    virtual void Dump() {}
    double const GetDelta() {return delta_;}
    int const GetNMembers() {return n_members_;}
    int const GetNPosit() {return n_posit_;}
    bool const GetPositFlag() {return posit_flag_;}
    void SetPositFlag(bool p_flag) {posit_flag_ = p_flag;}
    void SetNPosit(int n_pos) {n_posit_ = n_pos;}
    virtual int GetCount() {return 0;}
    virtual void Configurator() {}
    virtual void WriteOutputs(std::string run_name) {}
    virtual void WritePosits() {}
    virtual void ReadPosits() {}
    virtual void InitOutputFile(std::string run_name);
    virtual void InitInputFile(std::string in_file, std::ios::streampos beg);
    virtual int OutputIsOpen(){ return oposit_file_.is_open(); }
    virtual int InputIsOpen(){ return iposit_file_.is_open(); }
    virtual void ClosePosit(); 
};

template <typename T>
class Species : public SpeciesBase {
  protected:
    std::vector<T*> members_;
  public:
    Species() {}
    // Initialize function for setting it up on the first pass
    virtual void InitConfig(system_parameters *params, space_struct *space, long seed) {
      SpeciesBase::InitConfig(params, space, seed);
    }
    // Configurator function must be overridden
    virtual void Configurator() {
      error_exit("ERROR, species needs to override configurator!\n");
    }
    Species(int n_members, system_parameters *params, space_struct *space, long seed) : SpeciesBase(n_members, params, space, seed) {}
    //Virtual functions
    virtual void Init();
    virtual void AddMember();
    virtual void AddMember(T* newmem);
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositions();
    virtual std::vector<Simple*> GetSimples();
    virtual double GetKineticEnergy();
    virtual double GetPotentialEnergy();
    virtual double GetTotalEnergy();
    virtual double GetDrMax();
    virtual void ZeroDr();
    virtual void ZeroForces();
    virtual void Dump();
    virtual int GetCount();
    virtual void WritePosits();
    virtual void ReadPosits();
    virtual void ScalePositions();
    virtual std::vector<T*>* GetMembers() {return &members_;}
};

template <typename T> 
void Species<T>::Init() {
  for (int i=0; i<n_members_; ++i) {
    T * member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
    member->Init();
    members_.push_back(member);
    delete member;
  }
}

template <typename T> 
void Species<T>::AddMember() {
  T* newmember = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
  newmember->Init();
  members_.push_back(newmember);
  n_members_++;
  delete newmember;
}

template <typename T> 
void Species<T>::AddMember(T* newmem) {
  members_.push_back(newmem);
  n_members_++;
}

template <typename T> 
void Species<T>::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    (*it)->Draw(graph_array);
}

template <typename T> 
void Species<T>::UpdatePositions() {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    (*it)->UpdatePosition();
}

template <typename T> 
std::vector<Simple*> Species<T>::GetSimples() {
  std::vector<Simple*> simples;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    std::vector<Simple*> sim_vec = (*it)->GetSimples();
    simples.insert(simples.end(), sim_vec.begin(), sim_vec.end());
  }
  return simples;
}

template <typename T> 
double Species<T>::GetKineticEnergy() {
  double ke=0;
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    ke+=(*it)->GetKineticEnergy();
  return ke;
}

template <typename T> 
double Species<T>::GetPotentialEnergy() {
  double pe=0;
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    pe+=(*it)->GetPotentialEnergy();
    /* The total potential energy is going to be half of the
     potential energy felt by each particle. Potential energy 
     is shared, so I need to avoid double counting. */
  return 0.5*pe;
}

template <typename T> 
double Species<T>::GetTotalEnergy() {
  double ke = GetKineticEnergy();
  double pe = GetPotentialEnergy();
  return ke+pe;
}

template <typename T> 
double Species<T>::GetDrMax() {
  double dr_mag2_max = 0.0;
  double dr_mag2 = 0.0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    dr_mag2 = (*it)->GetDrMax();
    if (dr_mag2 > dr_mag2_max)
      dr_mag2_max = dr_mag2;
  }
  return dr_mag2_max;
}

template <typename T> 
void Species<T>::ZeroDr() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    (*it)->ZeroDrTot();
  }
}

template <typename T> 
void Species<T>::ZeroForces() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    (*it)->ZeroForce();
  }
}

template <typename T> 
void Species<T>::Dump() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    (*it)->Dump();
  }
}

template <typename T> 
int Species<T>::GetCount() {
  int count = 0;
  for (auto it =members_.begin(); it!=members_.end(); ++it) {
    count += (*it)->GetCount();
  }
  return count;
}

template <typename T> 
void Species<T>::WritePosits() {
  int size = members_.size();
  oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
  for( auto& mem_it : members_)
    mem_it->WritePosit(oposit_file_);
}

template <typename T> 
void Species<T>::ReadPosits() {
  if (iposit_file_.eof()) return;
  int size;
  bool del = false;
  T *member;
  iposit_file_.read(reinterpret_cast<char*>(&size), sizeof(size));

  if (size != members_.size()) {
    member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
    members_.resize(size, member);
    del = true;
  }
  for( auto& mem_it : members_){
    mem_it->ReadPosit(iposit_file_);
  }
  if (del)
    delete member;
}

template <typename T> 
void Species<T>::ScalePositions() {
  for (auto mem_it : members_) {
    mem_it->ScalePosition();
  }
}

#endif // _SIMCORE_SPECIES_H_
