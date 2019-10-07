#ifndef _SIMCORE_PARSE_PARAMS_H_
#define _SIMCORE_PARSE_PARAMS_H_

#include "parameters.hpp"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <string>

system_parameters parse_system_params(YAML::Node &node) {
  system_parameters params;
  return params;
}

species_base_parameters *parse_species_params(YAML::const_iterator it) {
  return new filament_parameters;
}

#endif // _SIMCORE_PARSE_PARAMS_H_
