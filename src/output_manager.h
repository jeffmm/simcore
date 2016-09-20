#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

void grabber(int width, int height, char *fname, int framenum);

class OutputManager{
  private:
    int *i_step_;
    double tot_virial_[3][3];

    YAML::Node node_;
    system_parameters *params_;
    std::vector<graph_struct*> *graph_array_;
    std::vector<std::string> posit_files_;

    std::map<SID, SpeciesBase*> species_;
    std::map<SID, int> posit_spec_; //Map of SIDs that are going be read in and 
                                    //how often each species will be read in
    void InitPositInput();

    void CalcVirial();
    void WriteVirial();
    void WriteSpeciesPosits();

    void AddSpecie(SpeciesBase *spec);

  public:
    OutputManager();
    void Init(system_parameters *params, 
        std::vector<graph_struct*> *graph_array, int *i_step);
    void SetMovie(std::vector<std::string> posit_files) { posit_files_ = posit_files;}
    bool IsMovie() { return !posit_files_.empty(); }
    void WriteOutputs();
    void AddSpecies(std::vector<SpeciesBase*> *species);
    void ReadSpeciesPositions();
    void GetGraphicsStructure();
    void Close();
    void Clear();

};

#endif //_SIMCORE_OUTPUT_MANAGER_H_

