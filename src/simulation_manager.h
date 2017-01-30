#ifndef _SIMCORE_SIMULATION_MANAGER_H_
#define _SIMCORE_SIMULATION_MANAGER_H_

#include "simulation.h"
//#include "analysis_manager.h"
#include "auxiliary.h"
#include "yaml-cpp/yaml.h"

class SimulationManager {

  private:
    bool make_movie_ = false,
         run_analysis_ = false;
    unsigned int n_runs_ = 1, //Number of runs per parameters set
                 n_var_ = 1; //Number of parameter variations
    std::string default_param_file_ = "src/master_params.yaml", 
                run_name_ = "sc"; // simulation batch name
    std::vector<std::string> pfiles_;
    Simulation *sim_; //New sim created and destroyed for every set of parameters
    YAML::Node pnode_; // Main node to initialize pvector
    std::vector<YAML::Node> pvector_;
    system_parameters params_; 
    rng_properties rng_;
    void CheckAppendParams();
    void AppendParams(YAML::Node app_node);
    void LoadDefaultParams();
    void CountVariations();
    void CheckRandomParams();
    double GetRandomParam(std::string rtype, int min, int max);
    void GenerateParameters();
    void WriteParams();
    void RunSimulations();
    void ParseParams(std::string file_name);
    //AnalysisManager analyzer_;

  public:
    SimulationManager() { debug_trace = false; }
    void InitManager(std::string param_file);
    //void RunManager(std::vector<std::string> posit_files);
    void RunManager();
    void SetNRuns(int n_runs) { n_runs_ = n_runs; }
    void SetRunName(std::string run_name) { run_name_ = run_name; }
    //void RunMovieManager(std::vector<std::string> posit_files);
    //void RunAnalyses(std::vector<std::string> pfiles);
    //void RunAnalysis() { run_analysis_ = true; }
    //void MakeMovie() { make_movie_ = true; }
    void DebugMode() { debug_trace = true; }
};


#endif // _SIMCORE_SIMULATION_MANAGER_H_ 
