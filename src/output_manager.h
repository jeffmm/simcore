#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

void grabber(int width, int height, char *fname, int framenum);

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
    //void SetInputPositFile();
    void WriteOutputs();
    void WriteSpeciesPosits();
    void AddSpecie(SpeciesBase *spec);
    std::ofstream& GetPositFile(SID sid);

};


#endif //_SIMCORE_OUTPUT_MANAGER_H_

