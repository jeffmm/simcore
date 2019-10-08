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
  for (std::size_t i = 0; i < labels_.size(); ++i) {
    if (species_id::_from_string(labels_[i].first.c_str()) == +species_id::crosslink) {
      n_crosslinks_++;
    } else {
      n_species_++;
    }
  }
}

void ParamsParser::ParseSpeciesParameters() {
  for (auto it = sim_node_.begin(); it != sim_node_.end(); ++it) {
    if (it->second.IsMap()) {
      if (it->first.as<std::string>() == "species") {
        continue;
      }
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
  CheckDuplicateLabels();
}

system_parameters ParamsParser::ParseSystemParameters() {
  return parse_system_params(sim_node_);
}

species_base_parameters *
ParamsParser::GetNewSpeciesParameters(sid_label &slab) {
  for (auto it = species_node_.begin(); it != species_node_.end(); ++it) {
    if (it->first.as<std::string>() == slab.first) {
      if (it->second.IsMap() && slab.second == "species") {
        return parse_species_params(slab.first, slab.second, it, sim_node_);
      } else if (it->second.IsSequence()) {
        for (int i = 0; i < it->second.size(); ++i) {
          YAML::Node subnode = it->second[i];
          for (YAML::const_iterator jt = subnode.begin(); jt != subnode.end();
               ++jt) {
            if (jt->first.as<std::string>() == slab.second) {
              return parse_species_params(slab.first, slab.second, jt, sim_node_);
            }
          }
        }
      }
    }
  }
  Logger::Error("ParamsParser did not find species id/label combination %s/%s"
                " in species params YAML node!",
                slab.first.c_str(), slab.second.c_str());
  return nullptr;
}
