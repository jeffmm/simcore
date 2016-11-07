#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "object.h"
#include <iterator>
#include "potential_base.h"

class OutputManager;

class SpeciesBase {
  private:
    SID sid_;
  protected:
    int n_members_;
    bool kmc_update_;
    bool is_kmc_ = false;
    double delta_;
    system_parameters *params_;
    std::fstream oposit_file_;
    std::fstream iposit_file_;

    space_struct *space_;
    al_set *anchors_;
    rng_properties rng_;
    void SetSID(SID sid) {sid_=sid;}

    double direct_[3] = {0};  //Neumatic director
    double pol_direct_[3] = {0}; //Polar neumatic director
    double virial_[9] = {0};

  public:
    SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed) {
      sid_ = SID::none;
      n_members_ = n_members;
      params_ = params;
      space_ = space;
      is_kmc_ = false;
      kmc_update_ = false;
      rng_.init(seed);
      delta_ = params->delta;
      ClearThermo();
    }
    SpeciesBase() {
      sid_ = SID::none;
      is_kmc_ = false;
      kmc_update_ = false;
    }
    virtual ~SpeciesBase() {
      rng_.clear();
      //for (auto it=potentials_.begin(); it!=potentials_.end(); ++it)
        //delete it->second;
      //potentials_.clear();
    }
    SpeciesBase(const SpeciesBase& that) {
      sid_=that.sid_;
      n_members_ = that.n_members_;
      params_=that.params_;
      space_=that.space_;
      //potentials_=that.potentials_;
      is_kmc_ = that.is_kmc_;
      kmc_update_ = that.kmc_update_;
      rng_.init(gsl_rng_get(that.rng_.r));
      delta_=that.delta_;
    }
    SpeciesBase& operator=(SpeciesBase const& that) {
      sid_=that.sid_;
      n_members_ = that.n_members_;
      params_=that.params_;
      space_=that.space_;
      //potentials_=that.potentials_;
      is_kmc_=that.is_kmc_;
      kmc_update_ = that.kmc_update_;
      rng_.init(gsl_rng_get(that.rng_.r));
      delta_=that.delta_;
      return *this;
    }
    virtual void UpdatePositions() {}
    virtual void UpdatePositionsMP() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void InitConfig(system_parameters *params, space_struct *space, al_set *pAnchors, long seed) {
      n_members_ = 0;
      params_ = params;
      space_ = space;
      anchors_ = pAnchors;
      //is_kmc_ = false;
      kmc_update_ = false;
      rng_.init(seed);
      delta_ = params->delta;
    }
    virtual void Init() {}
    virtual void ReInit(unsigned int const cid) {}
    //virtual void InitVirial() {}
    virtual double GetDrMax() {return 0.0;}
    virtual void ZeroDr() {}
    virtual void ZeroForces() {}
    virtual std::vector<Simple*> GetSimples() {
      std::vector<Simple*> sim;
      return sim;
    }
    virtual double GetKineticEnergy() {return 0;}
    virtual double GetPotentialEnergy() {return 0;}
    virtual double GetTotalEnergy() {return 0;}
    //virtual double const * const GetDirector(){ return direct_; }
    //virtual double const * const GetPolarDirector(){ return pol_direct_; }
    virtual void GetVirial(double *tot_virial){ 
      for (int i=0; i<params_->n_dim; i++)
        for (int j=i; j<params_->n_dim; j++){
          tot_virial[3*i+j] += virial_[3*i+j];
          tot_virial[3*j+i] = virial_[3*i+j];
        }
    }

    virtual void SetVirial(double* vl) {std::copy(vl, vl+9, virial_);}
    SID const GetSID() {return sid_;}
    bool IsKMC() {return is_kmc_;}
    bool const GetUpdate() {return kmc_update_;}
    void SetUpdate(bool update) {kmc_update_=update;}
    virtual void PrepKMC() {}
    virtual void StepKMC() {}
    virtual void DumpKMC() {}
    virtual void Dump() {}
    double const GetDelta() {return delta_;}
    int const GetNMembers() {return n_members_;}
    //std::vector<potential_pair> GetPotentials() {return potentials_;}
    virtual int GetCount() {return 0;}
    //std::vector<potential_pair> GetPotentials() {return potentials_;}
    virtual void Configurator() {}
    virtual void WriteOutputs(std::string run_name) {}
    virtual void WritePosits() {}
    virtual void ReadPosits() {}
    virtual void InitOutputFile(std::string run_name) {
      std::string sid_str = SIDToString(sid_);
      std::string file_name = run_name + "_" + sid_str + ".posit";
      std::cout<<"Posit file name is "<< run_name <<"_"<< sid_str <<" \n";
      oposit_file_.open(file_name, std::ios::out | std::ios::binary ); 
      if (!oposit_file_.is_open())
        std::cout<<"Output "<< file_name <<" file did not open\n";
      else{
        int size = sid_str.size();
        oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(int));
        oposit_file_.write(sid_str.c_str(), sid_str.size());
        oposit_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
        oposit_file_.write(reinterpret_cast<char*> (&params_->n_posit), sizeof(int));
      }
    }

    virtual void InitInputFile(std::string in_file, std::ios::streampos beg){
      iposit_file_.open(in_file, std::ios::binary | std::ios::in );
      iposit_file_.seekg(beg);
    }

    virtual int OutputIsOpen(){ return oposit_file_.is_open(); }
    virtual int InputIsOpen(){ return iposit_file_.is_open(); }
    virtual void Close(){ 
      if (oposit_file_.is_open())
        oposit_file_.close(); 
      if (iposit_file_.is_open())
        iposit_file_.close(); 
    }
    virtual std::vector<std::pair<unsigned int, unsigned int>> GetInternalPairs() {
      std::vector<std::pair<unsigned int, unsigned int>> retval;
      return retval;
    }
    virtual void ClearThermo(){}
    virtual double const * const GetDirector(){return nullptr;}
    virtual double const * const GetPolarDirector(){return nullptr;}
};

template <typename T>
class Species : public SpeciesBase {
  protected:
    std::vector<T*> members_;
  public:
    // Default constructor, needed for species factories
    Species() {}

    // Initialize function for setting it up on the first pass
    virtual void InitConfig(system_parameters *params, space_struct *space, al_set* pAnchors, long seed) {
      SpeciesBase::InitConfig(params, space, pAnchors, seed);
    }

    // Configurator function must be overridden
    virtual void Configurator() {
      printf("ERROR, needs to override!\n");
      exit(1);
    }

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
    //Virtual functions
    virtual void Init() {
      for (int i=0; i<n_members_; ++i) {
        T * member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
        member->Init();
        members_.push_back(member);
        delete member;
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

    virtual void AddMember() {
      T* newmember = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
      newmember->Init();
      members_.push_back(newmember);
      n_members_++;
      delete newmember;
    }

    virtual void AddMember(T* newmem) {
      members_.push_back(newmem);
      n_members_++;
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
      double dr_mag2 = 0.0;
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        dr_mag2 = (*it)->GetDrMax();
        if (dr_mag2 > dr_mag2_max)
          dr_mag2_max = dr_mag2;
      }
      //printf("dr2_max: %2.6f\n",dr_mag2_max);
      //if (dr_mag2_max > 6000) error_exit("!\n");
      return dr_mag2_max;
    }
    virtual void ZeroDr() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->ZeroDrTot();
      }
    }
    virtual void ZeroForces() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->ZeroForce();
      }
      ClearThermo();
    }
    virtual void Dump() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->Dump();
      }
    }
    virtual int GetCount() {
      int count = 0;
      for (auto it =members_.begin(); it!=members_.end(); ++it) {
        count += (*it)->GetCount();
      }
      return count;
    }
    
    virtual void WritePosits() {
      int size = members_.size();
      oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(size));
      for( auto& mem_it : members_)
        mem_it->WritePosit(oposit_file_);
    }

    virtual void ReadPosits() {
      T * member = new T(params_, space_, gsl_rng_get(rng_.r), GetSID());
      int size;
      iposit_file_.read(reinterpret_cast<char*>(&size), sizeof(size));
      members_.resize(size, member);
      for( auto& mem_it : members_){
        mem_it->ReadPosit(iposit_file_);
      }
      delete member;
    }

    virtual std::vector<std::pair<SID, SID>> GetInternalInteractionSIDs() {
      std::vector<std::pair<SID, SID>> sid_pairs;

      return sid_pairs;
    }
    virtual std::vector<std::pair<unsigned int, unsigned int>> GetInternalPairs() {
      std::vector<std::pair<unsigned int, unsigned int>> retval;

      for (auto it = members_.begin(); it != members_.end(); ++it) {
        auto mem_internal_pairs = (*it)->GetInternalPairs();
        retval.insert(retval.end(), mem_internal_pairs.begin(), mem_internal_pairs.end());
      }

      return retval;
    }

    virtual std::vector<T*>* GetMembers() {return &members_;}

    //virtual void GetVirial(double *tot_virial){ 
      //for (int i=0; i<params_->n_dim; i++)
        //for (int j=i; j<params_->n_dim; j++)
          //tot_virial[3*i+j] = tot_virial[3*j+i] += virial_[3*i+j];
    //}

    virtual void ClearThermo(){
      std::fill(direct_, direct_+3, 0);
      std::fill(pol_direct_, pol_direct_+3, 0);
      std::fill(virial_, virial_+9, 0);
      //memset(virial_, 0, sizeof(virial_));
    }
    virtual double const * const GetDirector(){
      int n_dim = params_->n_dim;
      double const *u;
      double north[] ={0,0,0};
      //Pick direction for upper hemisphere based on dimensionality
      //This is done because the director is a psuedo vector
      (n_dim == 2) ? (north[2] = 1) : (north[3] = 1);

      for (auto mem_it: members_){
        u = mem_it->GetOrientation();
        for (int i=0; i<n_dim; i++)
          direct_[i] += SIGN( u[i], dot_product(n_dim, north, u));
      }
      normalize_vector(direct_, n_dim);
      return direct_;
    }

    virtual double const * const GetPolarDirector(){
      int n_dim = params_->n_dim;
      double const *u;

      for (auto mem_it: members_){
        u = mem_it->GetOrientation();
        for (int i=0; i<n_dim; i++)
          pol_direct_[i] += u[i];
      }
      for (int i=0; i<n_dim; i++)
        pol_direct_[i] /= members_.size();
      //normalize_vector(pol_direct_, n_dim);
      return pol_direct_;
    }

    //virtual void InitVirial() {
    //  //memset(virial_, 0, sizeof(double)*9);
    //  std::fill(virial_, virial_+9, 0);
    //  for (auto mem : members_)
    //    mem->InitVirial(virial_);
    //}
};

#endif // _SIMCORE_SPECIES_H_
