#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

class OutputManager{
  private:
    int *i_step_,
        n_posit_;
    std::string run_name_;
    system_parameters *params_;
    std::vector<std::string> posit_files_;
    std::vector<SpeciesBase*> *species_;
    void InitPositInput();
    void WritePosits();

  public:
    OutputManager();
    ~OutputManager() {}
    void Init(system_parameters *params, 
              std::vector<SpeciesBase*> *species,
              int *i_step, std::string run_name);
    void SetMovie(std::vector<std::string> posit_files) {
      posit_files_ = posit_files; }
    void WriteOutputs();
    void ReadPosits();
    void Close();
};

#endif //_SIMCORE_OUTPUT_MANAGER_H_

