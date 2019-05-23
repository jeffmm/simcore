#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"
#include "yaml-cpp/yaml.h"

class SpeciesBase {
  private:
    species_id sid_;
  protected:
    int n_members_ = 0,
        spec_file_iterator_ = -1;
    std::string spec_name_;
    system_parameters *params_;
    species_parameters *sparams_;
    space_struct *space_;
    RNG rng_;
    std::fstream oposit_file_;
    std::fstream iposit_file_;
    std::fstream ospec_file_;
    std::fstream ispec_file_;
    std::string checkpoint_file_;
    std::vector<std::string> spec_file_names_;

  public:
    SpeciesBase() {}
    SpeciesBase(system_parameters *params, space_struct *space, long seed) {Init(params,space,seed);}
    void SetSID(species_id sid) {sid_=sid; spec_name_ = sid_._to_string();}
    virtual void UpdatePositions() {}
    virtual void Draw(std::vector<graph_struct*> * graph_array) {}
    virtual void Init(system_parameters *params, space_struct *space, long seed);
    virtual void ZeroForces() {}
    virtual std::vector<Object*> GetInteractors();
    virtual std::vector<Object*> GetLastInteractors();
    virtual double GetPotentialEnergy() {return 0;}
    virtual void ScalePositions() {}
    virtual void AddMember() {}
    virtual void SetLastMemberPosition(double const * const pos) {}
    virtual void PopMember() {}
    virtual void PopAll() {}
    virtual double GetSpecLength() {}
    virtual double GetSpecDiameter() {}
    virtual void ArrangeMembers() {}
    virtual int CanOverlap() {return sparams_->overlap;}
    species_id const GetSID() {return sid_;}
    virtual void Report() {}
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
    virtual void ReadPositsFromSpecs() {}
    virtual void InitAnalysis() {}
    virtual void RunAnalysis() {}
    virtual void FinalizeAnalysis() {}
    virtual void InitOutputFiles(std::string run_name);
    virtual void InitPositFile(std::string run_name);
    virtual void InitSpecFile(std::string run_name);
    virtual void InitPositFileInput(std::string run_name);
    virtual void InitSpecFileInput(std::string run_name);
    virtual bool InitSpecFileInputFromFile(std::string run_name);
    virtual bool HandleEOF();
    virtual void InitInputFiles(std::string run_name, bool posits_only, bool with_reloads);
    virtual void InitCheckpoints(std::string run_name);
    virtual void LoadFromCheckpoints(std::string run_name,std::string checkpoint_run_name);
    virtual int OutputIsOpen(){ return oposit_file_.is_open(); }
    virtual int InputIsOpen(){ return iposit_file_.is_open(); }
    virtual void CloseFiles(); 
    virtual void CleanUp() {}
    virtual void Reserve() {}
    virtual double const GetVolume() {}
    virtual double const GetDrMax() {}
    virtual void ZeroDrTot() {}
    virtual void CustomInsert() {}
};

template<typename T>
class Species : public SpeciesBase {
  protected:
    std::vector<T> members_;
  public:
    Species() {}
    // Initialize function for setting it up on the first pass
    virtual void Init(system_parameters *params, space_struct *space, long seed) {
      SpeciesBase::Init(params, space, seed);
    }
    Species(system_parameters *params, space_struct *space, long seed) : SpeciesBase(params, space, seed) {}
    //Virtual functions
    virtual void AddMember();
    virtual void AddMember(T newmem);
    virtual void PopMember();
    virtual void PopAll();
    virtual void ArrangeMembers();
    virtual void CrystalArrangement();
    virtual void CenteredOrientedArrangement();
    void SetLastMemberPosition(double const * const pos);
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositions();
    virtual std::vector<Object*> GetInteractors();
    virtual std::vector<Object*> GetLastInteractors();
    virtual double GetPotentialEnergy();
    virtual void ZeroForces();
    virtual void Report();
    virtual int GetCount();
    virtual void WritePosits();
    virtual void WriteSpecs();
    virtual void WriteCheckpoints();
    virtual void ReadPosits();
    virtual void ReadPositsFromSpecs();
    virtual void ReadSpecs();
    virtual void ReadCheckpoints();
    virtual void ScalePositions();
    virtual void InitAnalysis() {}
    virtual void RunAnalysis() {}
    virtual void FinalizeAnalysis() {}
    virtual std::vector<T> * GetMembers() {return &members_;}
    virtual void CleanUp();
    virtual void Reserve();
    virtual double const GetVolume();
    virtual double const GetDrMax();
    virtual void ZeroDrTot();
    virtual double GetSpecLength();
    virtual double GetSpecDiameter();
    virtual void CustomInsert();
};

#include "species_templates.h"

#endif // _SIMCORE_SPECIES_H_
