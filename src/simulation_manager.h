#ifndef _CYTOSCORE_SIMULATION_MANAGER_H_
#define _CYTOSCORE_SIMULATION_MANAGER_H_

//#include "simulation.h"
#include <yaml-cpp/yaml.h>
#include "auxiliary.h"

class SimulationManager {

  private:
    int n_runs_,
        n_var_;
    std::string run_name_;
    std::string param_file_;
    //Simulation *sim_;
    system_parameters params_;
    rng_properties rng_;
    void ParseParams(std::string param_file);
    void GetSeed(std::string param_file);
    void CheckVariations();
    void RunSimulations();

  public:
    SimulationManager(); // Default to debug mode
    SimulationManager(std::string param_file, std::string run_name, int n_runs);

};



#endif // _CYTOSCORE_SIMULATION_MANAGER_H_ 
