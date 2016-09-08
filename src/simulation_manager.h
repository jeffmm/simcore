#ifndef _SIMCORE_SIMULATION_MANAGER_H_
#define _SIMCORE_SIMULATION_MANAGER_H_

#include "simulation.h"
//#include <yaml.h>
#include <yaml-cpp/yaml.h>
#include "auxiliary.h"

class SimulationManager {

  private:
    unsigned int n_runs_, //Number of runs per parameters set
                 n_var_; //Number of parameter variations
    std::string run_name_;
    std::string param_file_; //Name of param yaml file
    Simulation *sim_; //New sim created and destroyed for every set of parameters
    system_parameters *params_; //array of system parameter structures for each sim that will be run
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



#endif // _SIMCORE_SIMULATION_MANAGER_H_ 
