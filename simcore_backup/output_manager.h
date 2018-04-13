#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

class OutputManager{
  private:
    bool posit_flag_ =  false,
         spec_flag_ =  false,
         checkpoint_flag_ = false,
         thermo_flag_ = false,
         posits_only_ = false;
    int *i_step_,
        n_posit_,
        n_spec_,
        n_checkpoint_,
        n_thermo_;
    std::string run_name_;
    system_parameters *params_;
    std::vector<SpeciesBase*> *species_;
    space_struct *space_;
    void WritePosits();
    void WriteSpecs();
    void WriteCheckpoints();
    void WriteThermo();
    void InitThermo();
    std::fstream thermo_file_;

  public:
    OutputManager() {}
    void Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species,
              space_struct *space,
              int *i_step, std::string run_name,
              bool reading_inputs = false,
              bool posits_only = false);
    //void SetPosits(std::vector<std::string> posit_files) {
      //posit_files_ = posit_files; 
    //}
    int GetNPosit() {return n_posit_;}
    int GetNSpec() {return n_spec_;}
    int GetNCheckpoint() {return n_checkpoint_;}
    void WriteOutputs();
    void InitInputs();
    void ReadInputs();
    void Close();
};

#endif //_SIMCORE_OUTPUT_MANAGER_H_
