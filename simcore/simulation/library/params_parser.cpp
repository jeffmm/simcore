#include "params_parser.hpp"
#include "parse_params.hpp"

void ParamsParser::Init(YAML::Node sim_params) {
  sim_node_ = sim_params;
  ParseSpeciesParameters();
}
void ParamsParser::CheckDuplicateLabels() {
  for (std::size_t i = 0; i < labels_.size() - 1; ++i) {
    for (std::size_t j = i + 1; j < labels_.size(); ++j) {
      if (labels_[i] == labels_[j]) {
        Logger::Error("Duplicate species label found! %s %s",
                      labels_[i].first.c_str(), labels_[i].second.c_str());
      }
    }
  }
}

void ParamsParser::ParseSpeciesParameters() {
  for (auto it = sim_node_.begin(); it != sim_node_.end(); ++it) {
    if (it->second.IsMap()) {
      species_node_[it->first.as<std::string>()] = it->second;
      labels_.push_back(std::make_pair(it->first.as<std::string>(), "species"));
    } else if (it->second.IsSequence()) {
      species_node_[it->first.as<std::string>()] = it->second;
      for (int i = 0; i < it->second.size(); ++i) {
        YAML::Node subnode = it->second[i];
        for (YAML::const_iterator jt = subnode.begin(); jt != subnode.end();
             ++jt) {
          labels_.push_back(std::make_pair(it->first.as<std::string>(),
                                           jt->first.as<std::string>()));
        }
      }
    }
  }
  for (YAML::const_iterator it = species_node_.begin(); it != species_node_.end(); ++it) {
    species_params_.push_back(parse_species_params(it));
  }
}

system_parameters ParamsParser::ParseSystemParameters() {
  return parse_system_params(sim_node_);
}
