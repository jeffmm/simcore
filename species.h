#ifndef _SIMCORE_SPECIES_H_
#define _SIMCORE_SPECIES_H_

#include "auxiliary.h"
#include "object.h"

class Species {
  protected:
    int n_members_ = 0;
    std::string spec_name_;
    std::vector<Object*> members_;
    system_parameters *params_;
    space_struct *space_;
    species_parameters *sparams_;
    std::fstream oposit_file_;
    std::fstream iposit_file_;
    std::fstream ospec_file_;
    std::fstream ispec_file_;
    std::string checkpoint_file_;
    RNG rng_;

  public:
    Species(system_parameters *params, space_struct *space, long seed);
    void SetSpecName(std::string specname) {spec_name_ = specname;}
    virtual void Init(system_parameters *params, space_struct *space, long seed);
    virtual int CanOverlap() {return sparams_->overlap;}
    int const GetNMembers() {return n_members_;}
    int const GetNInsert() {return sparams_->num;}
    int const GetNPosit() {return sparams_->n_posit;}
    int const GetNSpec() {return sparams_->n_spec;}
    int const GetNCheckpoint() {return sparams_->n_checkpoint;}
    bool const GetPositFlag() {return sparams_->posit_flag;}
    bool const GetSpecFlag() {return sparams_->spec_flag;}
    bool const GetCheckpointFlag() {return sparams_->checkpoint_flag;}
    std::string GetInsertionType() {return sparams_->insertion_type;}
    virtual void WriteOutputs(std::string run_name) {}
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
    // ---------------
    virtual void AddMember();
    virtual void AddMember(Object* newmem);
    virtual void PopMember();
    virtual void PopAll();
    virtual void ArrangeMembers();
    virtual void CrystalArrangement();
    virtual void CenteredOrientedArrangement();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositions();
    virtual std::vector<Object*> GetInteractors();
    virtual double GetPotentialEnergy();
    virtual void ZeroForces();
    virtual void Report();
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
    virtual std::vector<Object*>* GetMembers() {return &members_;}
    virtual void CleanUp();
};

#endif // _SIMCORE_SPECIES_H_
