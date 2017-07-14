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
    virtual void PopAll() {}
    virtual void ArrangeMembers() {}
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
    std::string GetInsertionType() {return sparams_->insertion_type;}
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
    virtual void CleanUp() {}
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
    virtual void PopAll();
    virtual void ArrangeMembers();
    virtual void SimpleCrystalArrangement();
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
    virtual void CleanUp();
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
void Species<T>::PopAll() {
  while (n_members_ > 0) {
    PopMember();
  }
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
    return;
  }
  if (! iposit_file_.is_open()) {
    printf(" ERROR: Posit file unexpectedly not open! Exiting.\n");
    early_exit=true;
    return;
  }
  int size = -1;
  T *member;
  iposit_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  // Hacky workaround FIXME
  // This prevents strange errors that occasionally crop up when reading inputs
  if (size == -1) {
    early_exit = true;
    return;
  }
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
  int size = 0;
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
    printf("  EOF reached\n");
    early_exit = true;
    return;
  }
  if (! ispec_file_.is_open()) {
    printf(" ERROR: Spec file unexpectedly not open! Exiting.\n");
    early_exit = true;
    return;
  }
  int size = -1;
  T *member;
  ispec_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
  // Hacky workaround FIXME
  if (size == -1) {
    early_exit = true;
    return;
  }
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

template <typename T>
void Species<T>::CleanUp() {
  for (auto it=members_.begin(); it!=members_.end(); ++it)
    delete (*it);
}

template<typename T>
void Species<T>::ArrangeMembers() {
  if (GetInsertionType().compare("simple_crystal") == 0)
    SimpleCrystalArrangement();
  else
    warning("WARNING! Arrangement not recognized and ArrangeMembers not overwritten by species!\n");
}

template<typename T>
void Species<T>::SimpleCrystalArrangement() {
  double d = members_[0]->GetDiameter();
  double l = members_[0]->GetLength();
  int n_dim = params_->n_dim;
  double R = space_->radius;
  double pos[3] = {0,0,0};
  double u[3] = {0,0,0};
  u[n_dim-1] = 1;
  double nX_max = floor((2.0*R/d)-1);
  double nY_max = floor((2.0*R-d)/(l+d));
  double nZ_max = (n_dim == 3 ? nX_max : 1);
  double N = nX_max*nY_max*nZ_max;
  if (n_members_ > N) {
    error_exit("ERROR: Number of members in species exceeds maximum possible for crystal in system radius! Max possible: %d\n", (int) N);
  }
  double fraction = N/n_members_;
  if (fraction < 1) fraction = 1;
  if (l == 0) {
    if (n_dim == 2)
      nY_max = floor(pow(n_members_, 1.0/n_dim));
    else if (n_dim == 3)
      nY_max = ceil(pow(n_members_, 1.0/n_dim));
  }
  if (n_dim == 2) {
    nX_max = ceil(n_members_/nY_max);
  }
  else if (n_dim == 3) {
    nX_max = floor(sqrt(n_members_/nY_max));
    if (nX_max == 0) nX_max = 1;
    nZ_max = ceil(n_members_/(nX_max*nY_max));
  }
  for (int i=0; i<n_dim; ++i) {
    pos[i] = -R+0.5*d;
  }
  pos[n_dim-1] += 0.5*l;
  int inserted = 0;
  int insert_x = 0;
  int insert_y = 0;
  int insert_z = 0;
  bool shift = false;
  double diff_y = (2*R - nY_max*(l+d))/nY_max;
  double diff_x = (2*R - nX_max*d)/nX_max;
  double diff_z = (2*R - nZ_max*d)/nZ_max;
  int u0 = 1;
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    // Check crystal orientation type
    if (params_->uniform_crystal == 0) {
      // Random orientations
      u[n_dim-1] = (gsl_rng_uniform_int(rng_.r,2) == 0 ? 1 : -1);
    }
    else if (params_->uniform_crystal == 2) {
      // Stagger orientations
      u[n_dim-1] = -u[n_dim-1];
    }
    // Insert object
    (*it)->InsertAt(pos,u);
    inserted++;
    // Update next position
    pos[n_dim-1] += l+d+diff_y;
    if (++insert_y == nY_max) {
      if (params_->uniform_crystal == 2 && shift) {
        // Stagger orientations row-wise
        u[n_dim-1] = u0;
        u0 = -u0;
      }
      // Stagger positions column-wise
      shift = !shift;
      insert_y = 0;
      pos[n_dim-1] = -R+0.5*(l+d) + (shift ? 0.5*(l+d+diff_y) : 0);
      pos[0] += d + diff_x;
      if (++insert_x == nX_max) {
        if (inserted < n_members_ && n_dim == 2) {
          error_exit("ERROR! Ran out of room while arranging crystal in 2D! Arranged %d/%d\n",inserted,n_members_);
        }
        else if (n_dim == 2) continue;
        insert_x = 0;
        pos[0] = -R+0.5*d;
        pos[1] += d + diff_z;
        if (++insert_z == nZ_max && inserted < n_members_) {
          error_exit("ERROR! Ran out of room while arranging crystal in 3D! Arranged %d/%d\n",inserted,n_members_);
        }
      }
    }
  }
}

#endif // _SIMCORE_SPECIES_H_
