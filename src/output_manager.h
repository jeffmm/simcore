#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"

void grabber(int width, int height, char *fname, int framenum);

typedef YAML::iterator yaml_it;

class OutputManager{
  private:
    int *i_step_;

    double tot_energy_;
    double tot_virial_[9];
    //std::vector< std::vector<double> > tot_virial_;
    double tot_direct_[3]; 
    double tot_pol_direct_[3];

    std::string run_name_;
    std::fstream thermo_file_;

    YAML::Node node_;
    system_parameters *params_;
    std::vector<graph_struct*> *graph_array_;
    std::vector<std::string> posit_files_;

    std::map<SID, SpeciesBase*> species_;
    std::map<SID, int> posit_spec_; //Map of SIDs that are going be read in and 
                                    //how often each species will be read in
    void InitPositInput();

    void CalcOutputs();
    void CalcSpeciesOutputs( SID sid, YAML::Node *data_type);


    void WriteSpeciesPosits();

    void WriteTotalThermo();
    void WriteSpeciesThermo();

    void AddSpecie(SpeciesBase *spec);

  public:
    OutputManager();
    ~OutputManager() {}
    void Init(system_parameters *params, 
        std::vector<graph_struct*> *graph_array, 
        int *i_step, std::string run_name);
    void MakeHeaders();
    void SetMovie(std::vector<std::string> posit_files) {
      posit_files_ = posit_files; }
    bool IsMovie() { return !posit_files_.empty(); }
    void WriteOutputs();
    void AddSpecies(std::vector<SpeciesBase*> *species);
    void ReadSpeciesPositions();
    void GetGraphicsStructure();
    void Close();
    void Clear();
};

#endif //_SIMCORE_OUTPUT_MANAGER_H_

