#include "simulation_manager.hpp"

bool early_exit;

/****************************************
   ::InitManaager::
   Initialize SimulationManager RNG and variables
   *************************************/
void SimulationManager::InitManager(run_options run_opts) {
  run_opts_ = run_opts;
  pnode_ = YAML::LoadFile(run_opts_.param_file);

  // Load parameter file specified run options and flags
  if (pnode_["n_runs"] && pnode_["n_runs"].size() == 0) {
    n_runs_ = pnode_["n_runs"].as<int>();
  }
  if (pnode_["n_random"] && pnode_["n_random"].size() == 0) {
    n_random_ = pnode_["n_random"].as<int>();
  }
  if (pnode_["run_name"] && pnode_["run_name"].size() == 0) {
    run_name_ = pnode_["run_name"].as<std::string>();
  }
  // Initialize logging
  InitLogger();
  long seed = 7143961348914;
  if (pnode_["seed"] && pnode_["seed"].size() == 0) {
    seed = pnode_["seed"].as<long>();
  } else {
    Logger::Warning("Default seed not overwritten!");
  }

  // Check command-line flags and run options.
  // Prefer command-line options over param file values.
  if (run_opts_.n_run_flag) {
    n_runs_ = run_opts_.n_runs;
  }
  if (run_opts_.run_name_flag) {
    run_name_ = run_opts_.run_name;
  }
  if (run_opts_.graphics_flag) {
    pnode_["graph_flag"] = 1;
    pnode_["n_graph"] = run_opts_.n_graph;
  }
  if ((run_opts_.make_movie || run_opts_.analysis_flag ||
       run_opts_.reduce_flag) &&
      n_runs_ > 1) {
    Logger::Error("Attempted to run movies/analysis on multiple files.");
    // Logger::Error("Attempted to run movies/analysis on multiple files.");
  }
  if (run_opts_.auto_graph) {
    pnode_["auto_graph"] = 1;
  }
  RNG::SetSeed(seed);
  rng_ = new RNG;

  // Check for appendable parameter files
  CheckAppendParams();
  // Load default parameters
  LoadDefaultParams();
  // Check for parameters that require randomized sequences for polynomial chaos
  // theory
  CheckRandomParams();
  // Count number of variations from parameter sequences
  CountVariations();
  // Generate n_var_ parameter nodes of unique parameter combinations
  GenerateParameters();
}

/* Initialize output logging. Requires that simulation run_name_ be initialized
 */
void SimulationManager::InitLogger() {
  std::ostringstream log_name;
  log_name << run_name_;
  if (run_opts_.load_checkpoint) {
    std::ostringstream nload;
    if (pnode_["n_load"]) {
      nload << std::setw(3) << std::setfill('0')
            << pnode_["n_load"].as<int>() + 1;
    } else {
      nload << std::setw(3) << std::setfill('0') << 1;
    }
    if (log_name.str().find("reload") == std::string::npos) {
      log_name << "_reload" << nload.str();
    } else {
      size_t pos = log_name.tellp();
      log_name.seekp(pos - 3);
      log_name << nload.str();
    }
  }
  Logger::SetOutput((log_name.str() + ".log").c_str());
}

/****************************************
   ::RunManaager::
   Main control sequence for SimulationManager
   *************************************/
void SimulationManager::RunManager() {
  // Write parameters to individual files
  WriteParams();

  if (run_opts_.blank_flag) {
    // Blank run -- we only write the parameter files without running a
    // simulation.
    return;
  } else if (run_opts_.analysis_flag || run_opts_.make_movie ||
             run_opts_.reduce_flag) {
    // Process the output files associated with the param file
    ProcessOutputs();
  } else {
    // Run simulations from parameter files
    RunSimulations();
  }
}

/****************************************
   ::CheckAppendParams::
   Check for existing "append_file" parameter in the main parameter node and if
   so treat values as file names to addition parameter sets and append
   corresponding parameter sets to the main parameter node.
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
      for (YAML::const_iterator it = app_params.begin(); it != app_params.end();
           ++it) {
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
  for (YAML::const_iterator it = app_node.begin(); it != app_node.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (!pnode_[param_name]) {
      if (it->second.IsMap()) {
        for (YAML::const_iterator jt = it->second.begin();
             jt != it->second.end(); ++jt) {
          pnode_[param_name][jt->first.as<std::string>()] = jt->second;
        }
      } else {
        pnode_[param_name] = it->second;
      }
    }
  }
}

/****************************************
   ::LoadDefaultParams::
   Open new YAML::Node from master parameter file (same one used in
   simcore_config) and load defaults into main parameter node if that node has
   not already initialized a value.
   *************************************/
void SimulationManager::LoadDefaultParams() {
  // YAML::Node defaults = YAML::LoadFile(params_.default_param_file);
#include "default_params.hpp"
  for (YAML::const_iterator it = default_config.begin();
       it != default_config.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (!pnode_[param_name]) {
      pnode_[param_name] = default_config[param_name];
    } else if (it->second.IsMap()) {
      for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
           ++jt) {
        std::string sub_param = jt->first.as<std::string>();
        if (!pnode_[param_name][sub_param]) {
          pnode_[param_name][sub_param] = default_config[param_name][sub_param];
        }
      }
    }
  }
}

/****************************************
   ::CheckRandomParam::
   Generates a random number N (either real or int) in the range (min,max)
   exclusive and returns N or 10^N, all depending on the random parameter
   generation type rtype, which can be R, RINT, or RLOG.
   *************************************/
void SimulationManager::CheckRandomParams() {
  std::string rtype;
  int n_params;
  double min, max;
  for (YAML::const_iterator it = pnode_.begin(); it != pnode_.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (it->second.IsSequence() && it->second.size() == 3 &&
        (rtype = it->second[0].as<std::string>()).at(0) == 'R') {
      min = it->second[1].as<double>();
      max = it->second[2].as<double>();
      pnode_[it->first] = YAML::Load("[R]");
      for (int i = 0; i < n_random_; ++i) {
        pnode_[it->first].push_back(GetRandomParam(rtype, min, max));
      }
    } else if (it->second.IsMap()) {
      for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
           ++jt) {
        if (jt->second.IsSequence() && jt->second.size() == 3 &&
            (rtype = jt->second[0].as<std::string>()).at(0) == 'R') {
          min = jt->second[1].as<double>();
          max = jt->second[2].as<double>();
          pnode_[it->first][jt->first] = YAML::Load("[R]");
          for (int i = 0; i < n_random_; ++i) {
            pnode_[it->first][jt->first].push_back(
                GetRandomParam(rtype, min, max));
          }
        }
      }
    }
  }
}

/****************************************
   ::GetRandomParam::
   Generates a random number N (either real or int) in the range (min,max)
   exclusive and returns N or 10^N, all depending on the random parameter
   generation type rtype, which can be R, RINT, or RLOG.
   *************************************/
double SimulationManager::GetRandomParam(std::string rtype, double min,
                                         double max) {
  if (max == min) {
    Logger::Error(
        "Min and max value of parameter randomization sequence are equal.");
  }
  if (rtype.compare("R") == 0) {
    return (min + (max - min) * gsl_rng_uniform_pos(rng_->r));
  } else if (rtype.compare("RINT") == 0) {
    return (min + gsl_rng_uniform_int(rng_->r, max - min));
  } else if (rtype.compare("RLOG") == 0) {
    return pow(10.0, min + (max - min) * gsl_rng_uniform_pos(rng_->r));
  } else {
    Logger::Error("Parameter randomization type not recognized.");
  }
}

/****************************************
   ::CountVariations::
   Given a main parameter file with a parameter value that is a sequence, count
   the total possible combinations of unique single parameter values
   *************************************/
void SimulationManager::CountVariations() {
  for (YAML::const_iterator it = pnode_.begin(); it != pnode_.end(); ++it) {
    if (it->second.IsSequence() &&
        (it->second[0].as<std::string>()).at(0) != 'R') {
      n_var_ *= it->second.size();
    } else if (it->second.IsMap()) {
      for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
           ++jt) {
        if (jt->second.IsSequence() &&
            (jt->second[0].as<std::string>()).at(0) != 'R') {
          n_var_ *= jt->second.size();
        }
      }
    }
  }
  n_var_ *= n_random_;
  if ((run_opts_.make_movie || run_opts_.analysis_flag ||
       run_opts_.reduce_flag) &&
      n_var_ > 1) {
    Logger::Error("Attempted to run movies/analysis on multiple files.");
  }
  if (n_var_ > 1) {
    Logger::Info("Initializing batch %s of %d simulations with %d variations of"
                 " %d runs each",
                 run_name_.c_str(), n_var_ * n_runs_, n_var_, n_runs_);
  } else if (n_runs_ > 1) {
    Logger::Info("Initializing batch of %d runs of simulation %s", n_runs_,
                 run_name_.c_str());
  } else {
    Logger::Info("Initializing simulation %s", run_name_.c_str());
  }
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
  int i_var, j_var, k_var = n_var_;
  for (YAML::const_iterator it = pnode_.begin(); it != pnode_.end(); ++it) {
    if (it->second.IsSequence() &&
        (it->second[0].as<std::string>()).at(0) != 'R') {
      int s = it->second.size();
      k_var /= s;
      j_var = n_var_ / (k_var * s);
      i_var = 0;
      for (int j = 0; j < j_var; ++j) {
        for (int i_param = 0; i_param < s; ++i_param) {
          for (int k = 0; k < k_var; ++k) {
            pvector_[i_var++][it->first] = it->second[i_param];
          }
        }
      }
    } else if (it->second.IsMap()) {
      for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
           ++jt) {
        if (jt->second.IsSequence() &&
            (jt->second[0].as<std::string>()).at(0) != 'R') {
          int s = jt->second.size();
          k_var /= s;
          j_var = n_var_ / (k_var * s);
          i_var = 0;
          for (int j = 0; j < j_var; ++j) {
            for (int i_param = 0; i_param < s; ++i_param) {
              for (int k = 0; k < k_var; ++k) {
                pvector_[i_var++][it->first][jt->first] = jt->second[i_param];
              }
            }
          }
        } else {
          for (i_var = 0; i_var < n_var_; ++i_var) {
            pvector_[i_var][it->first][jt->first] = jt->second;
          }
        }
      }
    } else {
      for (i_var = 0; i_var < n_var_; ++i_var) {
        pvector_[i_var][it->first] = it->second;
      }
    }
  }

  // Now handle random parameters
  k_var = n_var_ / n_random_;
  for (YAML::const_iterator it = pnode_.begin(); it != pnode_.end(); ++it) {
    if (it->second.IsSequence() &&
        (it->second[0].as<std::string>()).at(0) == 'R') {
      i_var = 0;
      for (int k = 0; k < k_var; ++k) {
        for (int i = 0; i < n_random_; ++i) {
          pvector_[i_var++][it->first] = it->second[i + 1];
        }
      }
    } else if (it->second.IsMap()) {
      for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
           ++jt) {
        if (jt->second.IsSequence() &&
            (jt->second[0].as<std::string>()).at(0) == 'R') {
          i_var = 0;
          for (int k = 0; k < k_var; ++k) {
            for (int i = 0; i < n_random_; ++i) {
              pvector_[i_var++][it->first][jt->first] = jt->second[i + 1];
            }
          }
        }
      }
    }
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
  for (int i_var = 0; i_var < n_var_; ++i_var) {
    for (int i_run = 0; i_run < n_runs_; ++i_run) {
      /* Only write new seed if we have more than one run/variation as
         this will ensure that parameters generated using WriteParams
         can be rerun individually with the expected result, ie not
         generating a different seed than the one generated here */
      if (n_runs_ > 1 || n_var_ > 1) {
        pvector_[i_var]["seed"] = gsl_rng_get(rng_->r);
      }
      std::ostringstream var;
      std::ostringstream run;
      std::ostringstream file_name;
      // setting zero padding to 3 digits
      var << std::setw(3) << std::setfill('0') << i_var;
      run << std::setw(3) << std::setfill('0') << i_run;
      file_name << run_name_;
      // append variation and run numbers to run_name if greater than 1
      if (n_var_ > 1) {
        file_name << "_v" << var.str();
      }
      if (n_runs_ > 1) {
        file_name << "_r" << run.str();
        pvector_[i_var]["n_runs"] = 1;
      }
      if (n_random_ > 1) {
        pvector_[i_var]["n_random"] = 1;
      }
      if (run_opts_.reduce_flag) {
        std::string red_file_name = file_name.str() + "_reduced" +
                                    std::to_string(run_opts_.reduce_factor);
        pvector_[i_var]["run_name"] = red_file_name;
        if (pvector_[i_var]["reduced"] == 0) {
          pvector_[i_var]["reduced"] = run_opts_.reduce_factor;
        }
        red_file_name = red_file_name + "_params.yaml";
        std::ofstream pfile(red_file_name, std::ios_base::out);
        YAML::Emitter out;
        pfile << (out << pvector_[i_var]).c_str();
        pfile.close();
      }
      if (run_opts_.load_checkpoint) {
        pvector_[i_var]["load_checkpoint"] = 1;
        int n_load = pvector_[i_var]["n_load"].as<int>() + 1;
        pvector_[i_var]["n_load"] = n_load;
        pvector_[i_var]["checkpoint_run_name"] = file_name.str();
        if (file_name.str().find("reduce") != std::string::npos) {
          pvector_[i_var]["reload_reduce_switch"] = 1;
        }
        std::ostringstream nload;
        nload << std::setw(3) << std::setfill('0') << n_load;
        if (file_name.str().find("reload") == std::string::npos) {
          file_name << "_reload" << nload.str();
        } else {
          size_t pos = file_name.tellp();
          file_name.seekp(pos - 3);
          file_name << nload.str();
        }
      }
      pvector_[i_var]["run_name"] = file_name.str();
      file_name << "_params.yaml";
      pfiles_.push_back(file_name.str());
      std::ofstream pfile(file_name.str(), std::ios_base::out);
      YAML::Emitter out;
      pfile << (out << pvector_[i_var]).c_str();
      pfile.close();
    }
  }
}

/***************************************
   ::RunSimulations::
   Loop over the vector of parameter file strings and load it as
   a YAML::Node. Parse the parameters of that node using the
   parse_params function (that is generated automatically using
   simcore_config) and create (and delete) a new simulation using
   those parameters.
   *************************************/
void SimulationManager::RunSimulations() {
  int n_sims = n_var_ * n_runs_;
  int i_sim = 1;
  for (std::vector<std::string>::iterator it = pfiles_.begin();
       it != pfiles_.end(); ++it) {
    ParseParams(*it);
    sim_ = new Simulation;
    Logger::Info("Starting simulation: %s", params_.run_name.c_str());
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
#include "parse_params.hpp"
void SimulationManager::ParseParams(std::string file_name) {
  YAML::Node node = YAML::LoadFile(file_name);
  YAML::Emitter out;
  Logger::Info("Initializing simulation with parameters:\n%s",
               (out << node).c_str());
  parse_params(node, &params_);
  if (run_opts_.load_checkpoint) {
    params_.load_checkpoint = 1;
  }
}

void SimulationManager::ProcessOutputs() {
  ParseParams(pfiles_[0]);
  sim_ = new Simulation;
  sim_->ProcessOutputs(params_, run_opts_);
  delete sim_;
}
