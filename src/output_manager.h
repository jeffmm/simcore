#ifndef _OUTPUT_MANAGER_H_
#define _OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

class OutputManager{
  private:
    int *i_step_;
    system_parameters *params_;
    std::map<SID, std::ofstream> posit_files_;
    std::map<SID, SpeciesBase*> species_;
    //Vector of species pointers

  public:
    OutputManager();
    void Init(system_parameters *params, int *i_step);
    void WriteOutputs();
    void WriteSpeciesPosits();
    void AddSpecie(SpeciesBase *spec);
    std::ofstream& GetPositFile(SID sid);

};

#endif //_SIMCORE_BR_BEAD_H_

