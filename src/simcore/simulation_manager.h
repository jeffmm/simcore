#ifndef _SIMCORE_SIMULATION_MANAGER_H_
#define _SIMCORE_SIMULATION_MANAGER_H_

#include "simulation.h"
#include "yaml-cpp/yaml.h"

class SimulationManager {

  private:
    bool make_movie_ = false,
         run_analysis_ = false;
    int n_runs_ = 1, //Number of runs per parameters set
                 n_var_ = 1, //Number of parameter variations
                 n_random_ = 1; //Number of random params
    std::string run_name_ = "sc"; // simulation batch name
    std::vector<std::string> pfiles_;
    Simulation *sim_; //New sim created and destroyed for every set of parameters
    YAML::Node pnode_; // Main node to initialize pvector
    std::vector<YAML::Node> pvector_;
    system_parameters params_; 
    run_options run_opts_;
    RNG rng_;
    void CheckAppendParams();
    void AppendParams(YAML::Node app_node);
    void LoadDefaultParams();
    void CountVariations();
    void SetRandomParams();
    double GetRandomParam(std::string rtype, double min, double max);
    void CheckRandomParams();
    void GenerateParameters();
    void WriteParams();
    void RunSimulations();
    void ParseParams(std::string file_name);
    void ProcessOutputs();
    UNIT_TESTER;

  public:
    SimulationManager() { debug_trace = false; early_exit = false;}
    void InitManager(run_options run_opts);
    void RunManager();
    //template <class T> struct access_by;
    //template <class T> friend struct access_by;
};

#endif // _SIMCORE_SIMULATION_MANAGER_H_ 
