#include "simulation_manager.h"

bool debug_trace;

/****************************************
   ::InitManaager::
   Initialize SimulationManager RNG and variables
   *************************************/
void SimulationManager::InitManager(run_options run_opts) {
  run_opts_ = run_opts;
  pnode_ = YAML::LoadFile(run_opts_.param_file);
  if (pnode_["n_runs"] && pnode_["n_runs"].size()==0)
    n_runs_ = pnode_["n_runs"].as<int>();
  if (pnode_["n_random"] && pnode_["n_random"].size()==0)
    n_random_ = pnode_["n_random"].as<int>();
  if (pnode_["run_name"] && pnode_["run_name"].size()==0)
    run_name_ = pnode_["run_name"].as<std::string>();
  long seed = 7143961348914;
  if (pnode_["seed"] && pnode_["seed"].size() == 0)
    seed = pnode_["seed"].as<long>();
  else
    std::cout << "  WARNING: Default seed not overwritten!\n";

  // Prefer command-line options over param values for n_runs, run_name
  if (run_opts_.debug)
    debug_trace = true;
  if (run_opts_.n_flag)
    n_runs_ = run_opts_.n_runs;
  if (run_opts_.r_flag)
    run_name_ = run_opts_.run_name;
  rng_.init(seed);
}

/****************************************
   ::RunManaager::
   Main control sequence for SimulationManager
   *************************************/
void SimulationManager::RunManager() {
  // Check for appendable parameter files
  CheckAppendParams();
  // Load default parameters
  LoadDefaultParams();
  // Check for parameters that require randomized sequences for polynomial chaos theory
  CheckRandomParams();
  // Count number of variations from parameter sequences 
  CountVariations();
  // Generate n_var_ parameter nodes of unique parameter combinations
  GenerateParameters();
  // Write parameters to individual files
  WriteParams();
  // Run simulations from parameter files
  RunSimulations();
}

/****************************************
   ::CheckAppendParams::
   Check for existing "append_file" parameter in the main parameter node and if so
   treat values as file names to addition parameter sets and append corresponding
   parameter sets to the main parameter node.
   *************************************/
void SimulationManager::CheckAppendParams() {
  if (YAML::Node app_params = pnode_["append_file"]) {
    // Add single file
    if (!app_params.IsSequence()) {
      YAML::Node app = YAML::LoadFile(app_params.as<std::string>());
      AppendParams(app);
    }
    // Add multiple files
    else {
      for (YAML::const_iterator it=app_params.begin(); it!=app_params.end(); ++it) {
        YAML::Node app = YAML::LoadFile(it->as<std::string>());
        AppendParams(app);
      }
    }
    pnode_.remove("append_file");
  }
}

/****************************************
   ::AppendParams::
   Append parameters from file that generated YAML::Node app_node to the main
   parameter node pnode_
   *************************************/
void SimulationManager::AppendParams(YAML::Node app_node) {
  for (YAML::const_iterator it=app_node.begin(); it!= app_node.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (!pnode_[param_name])
      if (it->second.IsMap())
        for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
          pnode_[param_name][jt->first.as<std::string>()] = jt->second;
      else
        pnode_[param_name] = it->second;
  }
}

/****************************************
   ::LoadDefaultParams::
   Open new YAML::Node from master parameter file (same one used in simcore_config)
   and load defaults into main parameter node if that node has not already
   initialized a value.
   *************************************/
void SimulationManager::LoadDefaultParams() {
  YAML::Node defaults = YAML::LoadFile(default_param_file_);
  for (YAML::const_iterator it=defaults.begin(); it!=defaults.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (!pnode_[param_name] && it->second.IsSequence())
        pnode_[param_name] = it->second[0];
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (!pnode_[param_name][jt->first.as<std::string>()] && jt->second.IsSequence())
          pnode_[param_name][jt->first.as<std::string>()] = jt->second[0];
  }
}

/****************************************
   ::CheckRandomParam::
   Generates a random number N (either real or int) in the range (min,max) exclusive
   and returns N or 10^N, all depending on the random parameter generation type
   rtype, which can be R, RINT, or RLOG.
   *************************************/
void SimulationManager::CheckRandomParams() {
  std::string rtype;
  int n_params;
  double min,max;
  for (YAML::const_iterator it=pnode_.begin(); it!=pnode_.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (it->second.IsSequence() && it->second.size()==3 && (rtype=it->second[0].as<std::string>()).at(0) == 'R') {
      min = it->second[1].as<double>();
      max = it->second[2].as<double>();
      pnode_[it->first] = YAML::Load("[R]");
      for (int i=0; i<n_random_; ++i)
        pnode_[it->first].push_back(GetRandomParam(rtype,min,max));
    }
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (jt->second.IsSequence() && jt->second.size()==3 && (rtype=jt->second[0].as<std::string>()).at(0) == 'R') {
          min = jt->second[1].as<double>();
          max = jt->second[2].as<double>();
          pnode_[it->first][jt->first] = YAML::Load("[R]");
          for (int i=0; i<n_random_; ++i)
            pnode_[it->first][jt->first].push_back(GetRandomParam(rtype,min,max));
        }
  }
}

/****************************************
   ::GetRandomParam::
   Generates a random number N (either real or int) in the range (min,max) exclusive
   and returns N or 10^N, all depending on the random parameter generation type
   rtype, which can be R, RINT, or RLOG.
   *************************************/
double SimulationManager::GetRandomParam(std::string rtype, int min, int max) {
  if (max == min) 
    error_exit("ERROR: Min and max value of parameter randomization sequence are equal.\n");
  if (rtype.compare("R") == 0) 
    return (min+(max-min)*gsl_rng_uniform_pos(rng_.r));
  else if (rtype.compare("RINT")==0) 
    return (min + gsl_rng_uniform_int(rng_.r,max-min));
  else if (rtype.compare("RLOG")==0) 
    return pow(10.0, min+(max-min)*gsl_rng_uniform_pos(rng_.r));
  else 
    error_exit("ERROR: Parameter randomization type not recognized.\n");
}

/****************************************
   ::CountVariations::
   Given a main parameter file with a parameter value that is a sequence, count
   the total possible combinations of unique single parameter values
   *************************************/
void SimulationManager::CountVariations() {
  for (YAML::const_iterator it=pnode_.begin(); it!=pnode_.end(); ++it) {
    if (it->second.IsSequence() && (it->second[0].as<std::string>()).at(0) != 'R')
      n_var_*=it->second.size();
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (jt->second.IsSequence() && (jt->second[0].as<std::string>()).at(0) != 'R')
          n_var_*=jt->second.size();
  }
  n_var_*=n_random_;
  if (n_var_ > 1)
    std::cout << "Initializing batch " << run_name_ << " of " << n_var_*n_runs_ << " simulations with " << n_var_ << " variations of " << n_runs_ << " runs each.\n";
  else if (n_runs_ > 1)
    std::cout << "Initializing batch of %d runs of simulation " << run_name_ <<".\n";
  else
    std::cout << "Initializing simulation " << run_name_ << ".\n";
}

/****************************************
   ::GenerateParameters::
   Generates n_var_ YAML::Nodes of parameter values each with a unique
   combination of single parameter values. This algorithm is of my own
   design and effectively disperses the parameter values over the range
   of parameter nodes that it shoud be applied in such a way as to make
   sure each parameter node has a unique combination, e.g. a parameter
   that is a sequence of 3 parameters will have each parameter
   dispersed over a third of the n_var_ parameter nodes, etc
    *************************************/
void SimulationManager::GenerateParameters() {
  pvector_.resize(n_var_);
  int i_var,j_var,k_var = n_var_;
  for (YAML::const_iterator it=pnode_.begin(); it!=pnode_.end(); ++it)
    if (it->second.IsSequence() && (it->second[0].as<std::string>()).at(0) != 'R') {
      unsigned int s = it->second.size();
      k_var /= s;
      j_var = n_var_ / (k_var * s) ;
      i_var = 0;
      for (int j=0; j<j_var; ++j)
        for (int i_param=0; i_param<s; ++i_param)
          for (int k=0; k<k_var; ++k)
            pvector_[i_var++][it->first] = it->second[i_param];
    }
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (jt->second.IsSequence() && (jt->second[0].as<std::string>()).at(0) != 'R') {
          unsigned int s = jt->second.size();
          k_var /= s;
          j_var = n_var_ / (k_var * s) ;
          i_var = 0;
          for (int j=0; j<j_var; ++j)
            for (int i_param=0; i_param<s; ++i_param)
              for (int k=0; k<k_var; ++k)
                pvector_[i_var++][it->first][jt->first] = jt->second[i_param];
        }
        else 
          for (i_var=0; i_var<n_var_; ++i_var) 
            pvector_[i_var][it->first][jt->first] = jt->second;
    else 
      for (i_var=0; i_var<n_var_; ++i_var)
        pvector_[i_var][it->first] = it->second;

  // Now handle random parameters
  k_var = n_var_/n_random_;
  for (YAML::const_iterator it=pnode_.begin(); it!=pnode_.end(); ++it)
    if (it->second.IsSequence() && (it->second[0].as<std::string>()).at(0) == 'R') {
      i_var=0;
      for (int k=0; k<k_var; ++k)
        for (int i=0; i<n_random_; ++i)
          pvector_[i_var++][it->first] = it->second[i+1];
    }
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (jt->second.IsSequence() && (jt->second[0].as<std::string>()).at(0) == 'R') {
          i_var=0;
          for (int k=0; k<k_var; ++k)
            for (int i=0; i<n_random_; ++i)
              pvector_[i_var++][it->first][jt->first] = jt->second[i+1];
        }
}

/****************************************
   ::WriteParams::
   Uses YAML::Emitter to write each of the parameter nodes to a YAML
   file, labeled by run name, variants, and runs, each with a unique
   simulation seed. Also stores each parameter file as a string in a
   vector to loop over in RunSimulations
   *************************************/
void SimulationManager::WriteParams() {
  for (int i_var=0; i_var<n_var_; ++i_var)
    for (int i_run=0; i_run<n_runs_; ++i_run) {
      /* Only write new seed if we have more than one run/variation as
         this will ensure that parameters generated using WriteParams
         can be rerun individually with the expected result, ie not
         generating a different seed than the one generated here */
      if (n_runs_ > 1 || n_var_ > 1)
        pvector_[i_var]["seed"]=gsl_rng_get(rng_.r);
      std::ostringstream var;
      std::ostringstream run;
      std::ostringstream file_name;
      // setting zero padding to 2 digits
      var << std::setw(2) << std::setfill('0') << i_var;
      run << std::setw(2) << std::setfill('0') << i_run;
      file_name << run_name_;
      // append variation and run numbers to run_name if greater than 1
      if (n_var_ > 1)
        file_name << "_v" << var.str();
      if (n_runs_ > 1) 
        file_name << "_r" << run.str();
      pvector_[i_var]["run_name"] = file_name.str();
      file_name << "_params.yaml";
      pfiles_.push_back(file_name.str());
      std::ofstream pfile(file_name.str(), std::ios_base::out);
      YAML::Emitter out;
      pfile << (out<<pvector_[i_var]).c_str();
      pfile.close();
    }
}

/****************************************
   ::RunSimulations::
   Loop over the vector of parameter file strings and load it as
   a YAML::Node. Parse the parameters of that node using the
   parse_params function (that is generated automatically using
   simcore_config) and create (and delete) a new simulation using
   those parameters. 
   *************************************/
void SimulationManager::RunSimulations() {
  int n_sims = n_var_*n_runs_;
  int i_sim = 1;
  for (std::vector<std::string>::iterator it=pfiles_.begin(); it!=pfiles_.end(); ++it) {
    ParseParams(*it);
    sim_ = new Simulation;
    std::cout << "\nRunning simulation " << params_.run_name<< " ("<<i_sim <<"/"<<n_sims<<")\n";
    sim_->Run(params_);
    delete sim_;
    i_sim++;
  }
}

/****************************************
   ::ParseParams::
   Load parameter file into a YAML::Node and parse the node for
   parameters in the system_parameters structure to initialize
   simulation parameters. Uses parse_params.h which is generated
   automatically using simcore_config.
   *************************************/
#include "parse_params.h"
void SimulationManager::ParseParams(std::string file_name) {
  YAML::Node node = YAML::LoadFile(file_name);
  YAML::Emitter out;
  std::cout << "Initializing simulation with parameters:\n" << (out<<node).c_str() << "\n";
  parse_params(node, &params_);
}

//void SimulationManager::RunMovieManager(std::vector<std::string> posit_files) {
  //Simulation *sim;
  //std::ostringstream title;

  ////FIXME will eventually have posit file read in this value
  //int i_var = 0; 
  //int i_run = 0;

  //InitVariations();
  //ParseParams(); //FIXME cannot take variations yet

  //title << run_name_;

  //if (n_var_ > 1)
    //title << "_v" << i_var+1;
  //if (n_runs_ > 1) 
    //title << "_r" << i_run+1;
  //params_[i_var].seed = gsl_rng_get(rng_.r);
  //PrintParams(params_[i_var], title.str());
  //sim = new Simulation;
  //sim->CreateMovie(params_[i_var], title.str(), posit_files);
  //delete sim;
  //title.str("");
  //title.clear();
//}

//void SimulationManager::RunAnalyses(std::vector<std::string> pfiles) {
  //// FIXME Decide how to use parameter variations with analysis arguments
  //// We'll just be using the first variation for now.
  //InitVariations();
  //ParseParams();
  //AnalysisManager aman;
  //std::ostringstream title;
  //// use first variation
  //int i_var = 0; 
  //// Set the run with a unique seed, in case we need random numbers in analysis.
  //params_[i_var].seed = gsl_rng_get(rng_.r);
  //analyzer_.RunAnalyses(params_[i_var], pfiles);
//}


