#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>
#include <ios>

class Configurator {
  private:
    YAML::Node node_;
    YAML::const_iterator subnode_;

    std::string master_ = "src/master_params.yaml",
                pfile_,
                pname_,
                pvalue_,
                ptype_,
                subname_,
                struct_name_,
                whitespace_;
    std::vector<std::string> pfiles_;

    std::ofstream parse_params_,
                  parameters_,
                  sub_parameters_;

    void WriteParameters();
    void WriteParam();
    void WriteParseParams();
    void WriteParseStructParams();
    void WriteParseSystemParams();
    void ParseParam(std::string iterator);

  public:
    Configurator() : 
        parse_params_("src/parse_params.h", std::ios_base::out),
        parameters_("src/parameters.h", std::ios_base::out) 
      {node_ = YAML::LoadFile(master_);}

    void ConfigureSimcore();
};

int main() {
  Configurator config;
  config.ConfigureSimcore();
  return 0;
}

void Configurator::ConfigureSimcore() {
  WriteParameters();
  WriteParseParams();
}

void Configurator::WriteParameters() {
  // Write header
  parameters_ << "#ifndef _SIMCORE_PARAMETERS_H_\n#define _SIMCORE_PARAMETERS_H_\n\nclass species_parameters;\n\n";
  // First parse sub parameters (species params, etc)
  std::vector<std::string> subparams;
  for (YAML::const_iterator it = node_.begin(); it!=node_.end(); ++it) {
    if (it->second.IsMap()) {
      std::ostringstream sub;
      if (it->first.as<std::string>().compare("species") == 0)
        parameters_ << "class " << it->first.as<std::string>() << "_parameters {\n  public:\n";
      else
        parameters_ << "class " << it->first.as<std::string>() << "_parameters : public species_parameters {\n  public:\n";
      for (subnode_ = it->second.begin(); subnode_!=it->second.end(); ++subnode_)
        WriteParam();
      parameters_ << "};\n\n";
      sub << it->first.as<std::string>() << "_parameters " << it->first.as<std::string>();
      subparams.push_back(sub.str());
    }
  }
  // Then parse remaining system parameters
  parameters_ << "class system_parameters {\n  public:\n";
  for (subnode_ = node_.begin(); subnode_!=node_.end(); ++subnode_)
    if (subnode_->second.IsSequence() && subnode_->second.size() == 2)
      WriteParam();
  // Add sub parameter structures
  for (std::vector<std::string>::iterator it=subparams.begin(); it!=subparams.end(); ++it)
    parameters_ << "    " << *it << ";\n";
  // Write end of file
  parameters_ << "};\n\n#endif // _SIMCORE_PARAMETERS_H_";
  parameters_.close();
}

void Configurator::WriteParam() {
  pname_ = subnode_->first.as<std::string>();
  pvalue_ = subnode_->second[0].as<std::string>();
  ptype_ = subnode_->second[1].as<std::string>();
  if (ptype_.compare("string") == 0) {
    parameters_ << "    std::string " << pname_ << " = \"" << pvalue_ << "\";\n";
  }
  else if (ptype_.compare("int") == 0 ||
      ptype_.compare("double") == 0 ||
      ptype_.compare("long") == 0 ||
      ptype_.compare("bool") == 0) {
    parameters_ << "    " << ptype_ << " " << pname_ << " = " << pvalue_ << ";\n";
  }
  else {
    std::cout << "  ERROR! Parameter type '" << ptype_ << "' not recognized.\n";
    exit(1);
  }
}

void Configurator::WriteParseParams() {

  // Write header
  parse_params_ << "#ifndef _SIMCORE_PARSE_PARAMS_H_\n#define _SIMCORE_PARSE_PARAMS_H_\n\n#include <iostream>\n#include <string>\n#include \"yaml-cpp/yaml.h\"\n#include \"parameters.h\"\n\nvoid parse_params(YAML::Node node, system_parameters *params) {\n\n  std::string param_name, param_value, struct_name;\n  for (YAML::const_iterator it=node.begin(); it!= node.end(); ++it) {\n    if (it->second.IsMap()) {\n      struct_name = it->first.as<std::string>();\n      if (false) {}\n";

  WriteParseStructParams();

  // Write transition to parse system params
  parse_params_ << "      else {\n        std::cout << \"  WARNING! Unrecognized struct parameter '\" << struct_name << \"'\\n\";\n      }\n    }\n    else {\n      param_name = it->first.as<std::string>();\n      if (false) {}\n";

  WriteParseSystemParams();

  // Write end of file
  parse_params_ << "      else {\n        std::cout << \"  WARNING! Unrecognized parameter '\" <<  param_name << \"'\\n\";\n      }\n    }\n  }\n}\n#endif // _SIMCORE_PARSE_PARAMS_H_";
  parse_params_.close();

}

void Configurator::WriteParseStructParams() {
  whitespace_ = "    ";
  for (YAML::const_iterator it=node_.begin(); it!=node_.end(); ++it) {
    if (it->second.IsMap()) {
      struct_name_ = it->first.as<std::string>();
      parse_params_ << "      else if (struct_name.compare(\""<< struct_name_ <<"\") == 0) {\n        for (YAML::const_iterator jt=it->second.begin(); jt!= it->second.end(); ++jt) {\n          param_name = jt->first.as<std::string>();\n          if (false) {}\n";

      struct_name_.append(".");
      for (subnode_ = it->second.begin(); subnode_!= it->second.end(); ++subnode_)
        ParseParam("jt");
      if (it->first.as<std::string>().compare("species")!=0)
        for (YAML::const_iterator jt=node_.begin(); jt!=node_.end(); ++jt)
          if (jt->second.IsMap() && jt->first.as<std::string>().compare("species")==0)
            for (subnode_ = jt->second.begin(); subnode_!=jt->second.end(); ++subnode_)
              ParseParam("jt");
      parse_params_ << "          else {\n            std::cout << \"  WARNING! Unrecognized \" << struct_name <<\" parameter: '\" << param_name << \"'\\n\";\n          }\n        }\n      }\n";
    }
  }
}

void Configurator::WriteParseSystemParams() {
  whitespace_ = "";
  struct_name_ = "";
  for (subnode_ = node_.begin(); subnode_!=node_.end(); ++subnode_)  {
    if (subnode_->second.IsSequence() && subnode_->second.size()==2) {
      ParseParam("it");
    }
  }
}

void Configurator::ParseParam(std::string iterator) {
  pname_ = subnode_->first.as<std::string>();
  ptype_ = subnode_->second[1].as<std::string>();
  if (ptype_.compare("string") == 0)
    ptype_ = "std::string";
  parse_params_ << whitespace_ << "      else if (param_name.compare(\"" << pname_ << "\")==0) {\n" << whitespace_ <<"        params->" << struct_name_ << pname_ << " = "<< iterator <<"->second.as<" << ptype_ << ">();\n" << whitespace_ <<"      }\n";
}

