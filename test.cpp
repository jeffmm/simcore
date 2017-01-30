#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "yaml-cpp/yaml.h"
#include "parameters.h"
#include "parse_params.h"

int main() {
  system_parameters parameters;
  YAML::Node params = YAML::LoadFile("params.yaml");
  std::string run_name_ = "sc";
  int n_var = 1;
  int n_runs = 1;
  if (params["n_runs"])
    n_runs = params["n_runs"].as<int>();
  printf ("num runs: %d\n",n_runs);

  // Check for appendable parameter files
  if (YAML::Node app_params = params["append_file"]) {
    if (!app_params.IsSequence()) {
      YAML::Node app = YAML::LoadFile(app_params.as<std::string>());
      for (YAML::const_iterator it=app.begin(); it!= app.end(); ++it) {
        std::string param_name = it->first.as<std::string>();
        if (!params[param_name])
          if (it->second.IsMap())
            for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
              params[param_name][jt->first.as<std::string>()] = jt->second;
          else
            params[param_name] = it->second;
      }
    }
    else {
      for (YAML::const_iterator it=app_params.begin(); it!=app_params.end(); ++it) {
        YAML::Node app = YAML::LoadFile(it->as<std::string>());
        for (YAML::const_iterator it=app.begin(); it!= app.end(); ++it) {
          std::string param_name = it->first.as<std::string>();
          if (!params[param_name])
            if (it->second.IsMap())
              for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
                params[param_name][jt->first.as<std::string>()] = jt->second;
            else
              params[param_name] = it->second;
        }
      }
    }
    params.remove("append_file");
  }

  // Load default parameters
  YAML::Node defaults = YAML::LoadFile("src/master_params.yaml");
  for (YAML::const_iterator it=defaults.begin(); it!=defaults.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    if (!params[param_name] && it->second.IsSequence())
        params[param_name] = it->second[0];
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (!params[param_name][jt->first.as<std::string>()] && jt->second.IsSequence())
          params[param_name][jt->first.as<std::string>()] = jt->second[0];
  }

  // Count number of variations from parameter sequences 
  for (YAML::const_iterator it=params.begin(); it!=params.end(); ++it) {
    if (it->second.IsSequence())
      n_var*=it->second.size();
    else if (it->second.IsMap())
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
        if (jt->second.IsSequence())
          n_var*=jt->second.size();
  }
  printf("num variations: %d\n",n_var);
  printf("num sims: %d\n", n_var*n_runs);

  // Create n_var parameter nodes of unique parameter combinations
  std::vector<YAML::Node> sub_params;
  sub_params.resize(n_var);
  int i_var,j_var,k_var = n_var;
  for (YAML::const_iterator it=params.begin(); it!=params.end(); ++it) {
    if (it->second.IsSequence()) {
      unsigned int s = it->second.size();
      k_var /= s;
      j_var = n_var / (k_var * s) ;
      i_var = 0;
      for (int j=0; j<j_var; ++j)
        for (int i_param=0; i_param<s; ++i_param)
          for (int k=0; k<k_var; ++k)
            sub_params[i_var++][it->first] = it->second[i_param];
    }
    else if (it->second.IsMap()) {
      for (YAML::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
        if (jt->second.IsSequence()) {
          unsigned int s = jt->second.size();
          k_var /= s;
          j_var = n_var / (k_var * s) ;
          i_var = 0;
          for (int j=0; j<j_var; ++j)
            for (int i_param=0; i_param<s; ++i_param)
              for (int k=0; k<k_var; ++k)
                sub_params[i_var++][it->first][jt->first] = jt->second[i_param];
        }
        else {
          for (i_var=0; i_var<n_var; ++i_var) 
            sub_params[i_var][it->first][jt->first] = jt->second;
        }
      }
    }
    else {
      for (i_var=0; i_var<n_var; ++i_var)
        sub_params[i_var][it->first] = it->second;
    }
  }

  // Output parameter files
  std::vector<std::string> files;
  for (i_var=0; i_var<n_var; ++i_var) {
    for (int i_run=0; i_run<n_runs; ++i_run) {
      sub_params[i_var]["seed"]=i_var;
      std::ostringstream var;
      std::ostringstream run;
      std::ostringstream file_name;
      // setting zero padding to 2 digits
      var << std::setw(2) << std::setfill('0') << i_var;
      run << std::setw(2) << std::setfill('0') << i_run;
      file_name << run_name_;
      if (n_var > 1)
        file_name << "_v" << var.str();
      if (n_runs > 1) 
        file_name << "_r" << run.str();
      file_name << "_params.yaml";
      files.push_back(file_name.str());
      std::ofstream file(file_name.str(), std::ios_base::out);
      YAML::Emitter out;
      file << (out<<sub_params[i_var]).c_str();
      file.close();
    }
  }

  for (std::vector<std::string>::iterator it=files.begin(); it!=files.end(); ++it) {
    YAML::Node node = YAML::LoadFile(*it);
    parse_params(node, &parameters);
    std::cout << "system radius: " << parameters.system_radius << "\n";
    std::cout << "filament driving: " << parameters.filament.driving_factor << "\n";
    std::cout << "filament insertion: " << parameters.filament.insertion_type << "\n";
  }

  return 0;
} 
