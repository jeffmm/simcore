#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

class OutputManager{
  private:
    bool posit_flag_,
         make_movie_ = false;
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
    void SetPosits(std::vector<std::string> posit_files) {
      make_movie_ = true;
      posit_files_ = posit_files; 
    }
    void WriteOutputs();
    bool IsMovie() {return make_movie_;}
    void ReadPosits();
    void Close();
};

#endif //_SIMCORE_OUTPUT_MANAGER_H_

