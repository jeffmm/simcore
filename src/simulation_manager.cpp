#include "simulation_manager.h"

SimulationManager::SimulationManager() {
  warning("Debug mode incomplete! Exiting!"); //FIXME
  exit(1);
}

SimulationManager::SimulationManager(std::string param_file, std::string run_name, int n_runs) {
  // Initialize parameters to default values
  params_.init(); 
  // Set run information
  YAML::Node node = YAML::LoadFile(param_file);
  YAML::iterator it=node.begin();
  std::cout << it->first.as<std::string>() << std::endl;
  std::cout << it->second[1].as<std::string>() << std::endl;
  std::cout << it->second.size() << std::endl;
  exit(0); //XXX
  n_runs_ = n_runs;
  param_file_ = param_file;
  run_name_ = run_name;
  // Check for variations of parameter values
  CheckVariations();
  // Initialize rng.
  GetSeed(param_file_);
  rng_.init(params_.seed);
  // Run all simulations
  RunSimulations();
}

void SimulationManager::GetSeed(std::string param_file) {
  YAML::Node node = YAML::LoadFile(param_file_);
  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if ( param_name.compare("seed") == 0 ) {
      std::string param_value = it->second.as<std::string>();
      params_.seed = atol(param_value.c_str()); 
    }
  }
}

void SimulationManager::CheckVariations() {
  std::cout << "Reading parameters from " << param_file_ << std::endl;
  YAML::Node node = YAML::LoadFile(param_file_);
  unsigned int n_var = 1;
  int n_arr = 0;
  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    if (it->second.size() > 1) {
      n_var *= it->second.size();
      n_arr++;
    }
  }
  if (n_var > 1) {
    std::cout << "Found " << n_var << " variations in " << param_file_ << ". Initializing " << n_var << " temporary parameter files of " << n_runs_ << " runs each for a total of " << n_runs_*n_var << " unique simulations." << std::endl;

    CreateVariations(n_arr);
  }
}

void SimulationManager::CreateVariations(int n_arr) {
  //int *size_arr = new int(n_arr);
  //int *i_arr = new int(n_arr);
  //int j_arr = 0;
  //for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    //if (it->second.size() > 1) {
      //size_arr[j_arr] = it->second.size();
      //i_arr[j_arr] = 0;
      //j_arr++;
    //}
  //}

  //int i_var = 0;
  //for (j_arr = 0; j_arr < n_arr; ++j_arr) {
    //for (int k_arr = 0; k_arr < size_arr[j_arr]; ++k_arr) {
      //std::ostringstream param_var_file;
      //param_var_file << param_file_ << i_var;
      //YAML::Node node = YAML::LoadFile(param_file_);
      //std::ofstream param_var((param_var_file.str()).c_str(), std::ios_base::out | std::ios_base::app);
      //for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
        //std::string param_name = it->first.as<std::string>();
        //std::string param_value = it->second.as<std::string>();
        //param_var << param_name << " : " << param_value << "\n";
      //}
    //}
  //}
}

void SimulationManager::ParseParams(std::string param_file) {
 
  std::cout << "Reading parameters from " << param_file << std::endl;
  YAML::Node node = YAML::LoadFile(param_file);

  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    std::string param_value = it->second.as<std::string>();
  
    // check for all parameters
#include "parse_params_body.cpp"
  }
  std::cout << std::endl;
}

// We should parallelize this in the future
// to submit all runs to available processors
void SimulationManager::RunSimulations() {
  for (int i_var=0; i_var<n_var_; ++i_var) {
    if (n_var_ > 1) {
      //ostringstream param_var = param_file << i_var;
      //ParseParams(param_var.str());
    }
    else {
      ParseParams(param_file_);
    }
  //sim_ = new Simulation(&params_);
  }
}
