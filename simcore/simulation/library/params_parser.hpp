#ifndef _SIMCORE_PARAMS_PARSER_H_
#define _SIMCORE_PARAMS_PARSER_H_

#include "auxiliary.hpp"
#include "yaml-cpp/yaml.h"

typedef std::pair<std::string, std::string> sid_label;

class ParamsParser {
private:
  YAML::Node sim_node_;
  YAML::Node species_node_;
  std::vector<sid_label> labels_;
  void CheckDuplicateLabels();
  void ParseSpeciesParameters();
  std::vector<species_parameters *> species_params_;

public:
  void Init(YAML::Node sim_params);
  system_parameters ParseSystemParameters();
  std::vector<species_parameters *> GetSpeciesParameters() {
    return species_params_;
  }
};

#endif
