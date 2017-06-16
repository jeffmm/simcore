#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"

class SpeciesBase {
  private:
    SID sid_;
  protected:
    int n_members_ = 0;
    system_parameters *params_;
    species_parameters *sparams_;
    std::fstream oposit_file_;
    std::fstream iposit_file_;
    std::fstream ospec_file_;
    std::fstream ispec_file_;
    std::string checkpoint_file_;
    space_struct *space_;
    rng_properties rng_;
    void SetSID(SID sid) {sid_=sid;}

  public:
    SpeciesBase() {sid_ = SID::none;}
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed);
    virtual void UpdatePositions() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init(system_parameters *params, space_struct *space, long seed);
    virtual double GetDrMax() {return 0.0;}
    virtual void ZeroDr() {}
    virtual void ZeroForces() {}
    virtual std::vector<Simple*> GetSimples();
    virtual double GetKineticEnergy() {return 0;}
    virtual double GetPotentialEnergy() {return 0;}
    virtual double GetTotalEnergy() {return 0;}
    virtual void ScalePositions() {}
    virtual void AddMember() {}
    virtual void PopMember() {}
    virtual int CanOverlap() {return sparams_->overlap;}
    SID const GetSID() {return sid_;}
    virtual void Dump() {}
    int const GetNMembers() {return n_members_;}
    int const GetNInsert() {return sparams_->num;}
    int const GetNPosit() {return sparams_->n_posit;}
    int const GetNSpec() {return sparams_->n_spec;}
    int const GetNCheckpoint() {return sparams_->n_checkpoint;}
    bool const GetPositFlag() {return sparams_->posit_flag;}
    bool const GetSpecFlag() {return sparams_->spec_flag;}
    bool const GetCheckpointFlag() {return sparams_->checkpoint_flag;}
    virtual int GetCount() {return 0;}
    virtual void WriteOutputs(std::string run_name) {}
    virtual void WritePosits() {}
    virtual void WriteSpecs() {}
    virtual void WriteCheckpoints() {}
    virtual void ReadSpecs() {}
    virtual void ReadCheckpoints() {}
    virtual void ReadPosits() {}
    virtual void InitAnalysis() {}
    virtual void RunAnalysis() {}
    virtual void FinalizeAnalysis() {}
    virtual void InitOutputFiles(std::string run_name);
    virtual void InitPositFile(std::string run_name);
    virtual void InitSpecFile(std::string run_name);
    virtual void InitPositFileInput(std::string run_name);
    virtual void InitSpecFileInput(std::string run_name);
    virtual void InitInputFiles(std::string run_name, bool posits_only);
    virtual void InitCheckpoints(std::string run_name);
    virtual int OutputIsOpen(){ return oposit_file_.is_open(); }
    virtual int InputIsOpen(){ return iposit_file_.is_open(); }
    virtual void CloseFiles(); 
};

template <typename T>
class Species : public SpeciesBase {
  protected:
    std::vector<T*> members_;
  public:
    Species() {}
    // Initialize function for setting it up on the first pass
    virtual void Init(system_parameters *params, space_struct *space, long seed) {
      SpeciesBase::Init(params, space, seed);
    }
    Species(int n_members, system_parameters *params, space_struct *space, long seed) : SpeciesBase(n_members, params, space, seed) {}
    //Virtual functions
    virtual void AddMember();
    virtual void AddMember(T* newmem);
    virtual void PopMember();
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
    virtual void WriteSpecs();
    virtual void WriteCheckpoints();
    virtual void ReadPosits();
    virtual void ReadSpecs();
    virtual void ReadCheckpoints();
    virtual void ScalePositions();
    virtual void InitAnalysis() {}
    virtual void RunAnalysis() {}
    virtual void FinalizeAnalysis() {}
    virtual std::vector<T*>* GetMembers() {return &members_;}
};

template <typename T> 
void Species<T>::AddMember() {
  T* newmember = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
  newmember->Init();
  //newmember->SetColor(sparams_->color, sparams_->draw_type);
  members_.push_back(newmember);
  n_members_++;
}

template <typename T> 
void Species<T>::AddMember(T* newmem) {
  members_.push_back(newmem);
  n_members_++;
}

template <typename T> 
void Species<T>::PopMember() {
  delete members_[n_members_-1];
  members_.pop_back();
  n_members_--;
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
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ZeroDrTot();
}

template <typename T> 
void Species<T>::ZeroForces() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ZeroForce();
}

template <typename T> 
void Species<T>::Dump() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->Dump();
}

template <typename T> 
int Species<T>::GetCount() {
  int count = 0;
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    count += (*it)->GetCount();
  return count;
}

template <typename T> 
void Species<T>::WritePosits() {
  int size = members_.size();
  oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WritePosit(oposit_file_);
}

template <typename T> 
void Species<T>::WriteSpecs() {
  int size = members_.size();
  ospec_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WriteSpec(ospec_file_);
}

template <typename T> 
void Species<T>::WriteCheckpoints() {
  int size = members_.size();
  std::fstream ocheck_file(checkpoint_file_,std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    std::cout<<"ERROR: Output "<< checkpoint_file_ <<" file did not open\n";
    exit(1);
  }
  ocheck_file.write(reinterpret_cast<char*>(&size), sizeof(size));
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->WriteCheckpoint(ocheck_file);
  ocheck_file.close();
}

template <typename T> 
void Species<T>::ReadPosits() {
  if (iposit_file_.eof()) {
    printf("  EOF reached\n");
    early_exit = true;
  }
  int size;
  T *member;
  iposit_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  if (size != n_members_) {
    member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
    members_.resize(size, member);
    delete member;
    n_members_ = size;
  }
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadPosit(iposit_file_);
}

template <typename T> 
void Species<T>::ReadCheckpoints() {
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    std::cout<<"  ERROR: Output "<< checkpoint_file_ <<" file did not open\n";
    exit(1);
  }
  int size;
  T *member;
  icheck_file.read(reinterpret_cast<char*>(&size), sizeof(size));

  member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
  members_.resize(size, member);
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadCheckpoint(icheck_file);
  icheck_file.close();
}

template <typename T> 
void Species<T>::ReadSpecs() {
  if (ispec_file_.eof()) {
    CloseFiles();
    printf("  EOF reached\n");
    early_exit = true;
  }
  int size;
  T *member;
  ispec_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  if (size != n_members_) {
    member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
    members_.resize(size, member);
    delete member;
    n_members_ = size;
  }
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ReadSpec(ispec_file_);
}

template <typename T> 
void Species<T>::ScalePositions() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) 
    (*it)->ScalePosition();
}

#endif // _SIMCORE_SPECIES_H_
