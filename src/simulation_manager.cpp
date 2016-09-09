#include "simulation_manager.h"

bool debug_trace;

SimulationManager::SimulationManager() {
  debug_trace = false;
}
SimulationManager::~SimulationManager() {
  // Delete system_parameters structure, initialized in InitVariations
  delete[] params_;
}

/****************************************
   ::RunManager::
   Runs (n_var_*n_runs_) simulations where n_var_ is the number of parameter
   variations given from combinatorics of lists of parameter values from the
   input parameter file (n_var_ = 1 if there are no lists of parameter values),
   and n_runs_ is the number of runs given at the command line by the -n flag
   or n_runs parameter value (1 by default)
   *************************************/
void SimulationManager::RunManager() {
  InitVariations();
  ParseParams();
  RunSimulations();
}

/****************************************
   ::DebugMode::
   TODO Create a series of simulation runs that captures validation data
   *************************************/
void SimulationManager::DebugMode() {
  //warning("Debug mode incomplete! Exiting!"); //FIXME
  //exit(1);
  debug_trace = true;
}


/****************************************
   ::InitManager::
   Initialize n_runs_ and run_name_ from parameter values and initialize RNG
   with seed from input file. 
   *************************************/
void SimulationManager::InitManager(std::string param_file) {
  n_runs_=1;  // default number of runs
  run_name_="sc"; // default run name
  param_file_ = param_file;
  YAML::Node node = YAML::LoadFile(param_file_);
  long seed = 71348958934175; // default seed... should be overwritten
  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    if (it->first.size()==0) {
      std::string param_name = it->first.as<std::string>();
      if ( param_name.compare("seed") == 0 ) {
        std::string param_value = it->second.as<std::string>();
        seed = atol(param_value.c_str()); 
      }
      else if (param_name.compare("n_runs") == 0) {
        std::string param_value = it->second.as<std::string>();
        n_runs_ = atoi(param_value.c_str());
      }
      else if (param_name.compare("run_name") == 0) {
        run_name_ = it->second.as<std::string>();
      }
    }
  }
  rng_.init(seed);
}

/****************************************
   ::InitVariations::
   Counts the number of variations in the parameter file, in case some
   parameters are given a list of values. This is simply a multiplication
   series of the size of each list of values, with single values having size
   one. Also initializes the system_parameters array with n_var_ params
   structures.
   *************************************/
void SimulationManager::InitVariations() {
  std::cout << "Reading parameters from " << param_file_ << std::endl;
  YAML::Node node = YAML::LoadFile(param_file_);
  n_var_ = 1;
  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    if (it->first.size() == 1) {
      for (YAML::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
        if (jt->second.size() > 1) {
          n_var_ *= jt->second.size();
        }
      }
    }
    else {
      if (it->second.size() > 1) {
        n_var_ *= it->second.size();
      }
    }
  }
  if (n_var_ > 1) {
    std::cout << "Found " << n_var_ << " variations in " << param_file_ << ". Initializing " << n_var_ << " temporary parameter files of " << n_runs_ << " runs each for a total of " << n_runs_*n_var_ << " unique simulations." << std::endl;
  }
  // Initialize each set of parameters with default values
  params_ = new system_parameters[n_var_];
}

/****************************************
   ::ParseParams::
   This function will parse the input parameter file and, if necessary,
   create n_var parameter files with single parameter values if the
   input parameter file included a list of values for any number of
   parameters. The algorithm is of my own design, but effectively
   creates all the combinatorics of unique parameter files from
   any lists of parameter values.
   *************************************/
void SimulationManager::ParseParams() {
  YAML::Node node = YAML::LoadFile(param_file_);
  unsigned int increment, i_var, j_var;
  unsigned int k_var = n_var_;
  std::string param_value, param_name;
  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    param_name = it->first.as<std::string>();
    if (it->second.size() > 1) {
      k_var /= it->second.size();
      j_var = n_var_/(k_var*it->second.size());
      i_var = 0;
      for (int j=0; j<j_var; ++j) {
        for (int i_param=0; i_param<it->second.size(); ++i_param) {
          for (int k=0; k<k_var; ++k) {
            param_value = it->second[i_param].as<std::string>();
            ParseParameter(param_name, param_value, i_var);
            i_var++;
          }
        }
      }
    }
    else {
      param_value = it->second.as<std::string>();
      for (i_var=0; i_var<n_var_; ++i_var)
        ParseParameter(param_name, param_value, i_var);
    }
  }
}

/****************************************
   ::ParseParameter::
   Checks single parameter value from YAML Node and initializes the parameter
   param_name with value for the i_var'th param structure
   *************************************/
void SimulationManager::ParseParameter(std::string param_name, std::string param_value, unsigned int i_var) {
#include "parse_params_body.h"
}

/****************************************
   ::PrintParams::
   Create new parameter file from params structure with prefix 'name'
   *************************************/
void SimulationManager::PrintParams(system_parameters params, std::string name) {
  std::ostringstream file_name;
  file_name << name << "-params.yaml";
  std::ofstream param_file((file_name.str()).c_str(), std::ios_base::out);
#include "print_params_body.h"
  param_file.close();
}

/****************************************
   ::RunSimulations::
   Run n_runs_ simulations for each unique parameter file, each with a unique
   seed ( though still determined from the SimulationManager RNG seed)
   TODO In the future we can do a simple parallelization to submit all runs to
   available processors
   *************************************/
void SimulationManager::RunSimulations() {
  Simulation *sim;
  std::ostringstream title;
  for (int i_var=0; i_var<n_var_; ++i_var) {
    for (int i_run=0; i_run<n_runs_; ++i_run) {
      // Set each run with a unique seed
      title << run_name_;
      if (n_var_ > 1)
        title << "-v" << i_var+1;
      if (n_runs_ > 1) 
        title << "-r" << i_run+1;
      params_[i_var].seed = gsl_rng_get(rng_.r);
      PrintParams(params_[i_var], title.str());
      sim = new Simulation;
      sim->Run(params_[i_var], title.str());
      delete sim;
      title.str("");
      title.clear();
    }
  }
}

void SimulationManager::SetNRuns(int n_runs) {
  n_runs_ = n_runs;
}

void SimulationManager::SetRunName(std::string run_name) {
  run_name_ = run_name;
}

void SimulationManager::RunMovieManager(std::string posit_file) {
  Simulation *sim;
  std::ostringstream title;

  //FIXME will eventually have posit file read in this value
  int i_var = 0; 
  int i_run = 0;

  InitVariations();
  ParseParams(); //FIXME cannot take variations yet

  title << run_name_;

  if (n_var_ > 1)
    title << "-v" << i_var+1;
  if (n_runs_ > 1) 
    title << "-r" << i_run+1;
  params_[i_var].seed = gsl_rng_get(rng_.r);
  PrintParams(params_[i_var], title.str());
  sim = new Simulation;
  sim->CreateMovie(params_[i_var], title.str(), posit_file);
  delete sim;
  title.str("");
  title.clear();
}
  

