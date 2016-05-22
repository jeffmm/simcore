#ifndef _CYTOSCORE_SIMULATION_MANAGER_H_
#define _CYTOSCORE_SIMULATION_MANAGER_H_

#include "simulation.h"
#include <yaml-cpp/yaml.h>
#include "auxiliary.h"

class SimulationManager {

  private:
    unsigned int n_runs_,
                 n_var_;
    std::string run_name_;
    std::string param_file_;
    Simulation *sim_;
    system_parameters *params_;
    rng_properties rng_;
    void ParseParams();
    void InitVariations();
    //void CreateVariations();
    void RunSimulations();
    void ParseParameter(std::string param_name, std::string param_value, unsigned int i_var);
    void PrintParams(system_parameters params, std::string name);

  public:
    SimulationManager();
    ~SimulationManager();
    void InitManager(std::string param_file);
    void DebugMode();
    void RunManager();
    void SetNRuns(int n_runs);
    void SetRunName(std::string run_name);
};



#endif // _CYTOSCORE_SIMULATION_MANAGER_H_ 
