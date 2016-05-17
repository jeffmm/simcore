#include <iostream>
#include <yaml-cpp/yaml.h>
#include <string.h>
#include <sstream>
#include <ios>
#include <fstream>

void init_files();
void write_params(std::string param_name, std::string param_value, std::string param_type);
void write_parse_params_body_cpp(std::string param_name, std::string param_type);
void write_parameters_h(std::string param_name, std::string param_value, std::string param_type);
std::string parse_params_snippet(std::string param_name, std::string param_type);
std::string parameters_snippet(std::string param_name, std::string param_value, std::string param_type);
void write_print_params_body_cpp(std::string param_name);
void cleanup_files();
void write_objects(std::string *object_list);

int main(int argc, char *argv[]) {

  std::string param_file, param_name, param_value, param_type;
  std::string *object_list;
  if (argc==1) 
    param_file = "src/params_master.yaml";
  else
    param_file = argv[1];

  init_files();
  YAML::Node node = YAML::LoadFile(param_file);
  for (YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    if (it->first.size()==0 && it->second.size() != 2) {
      std::cout << "  ERROR! Master Parameter file parameter values must be in the form of [ param_value, param_type ]\n";
      exit(1);
    }
    else if (it->first.size()==0) {
      param_name = it->first.as<std::string>();
      param_value = it->second[0].as<std::string>();
      param_type = it->second[1].as<std::string>();
      write_params(param_name, param_value, param_type);
    }
    else if (it->first.size()==1 && (it->first[0].as<std::string>()).compare("objects")==0) {
      object_list = new std::string[it->second.size()];
      for (unsigned int i=0; i<it->second.size(); ++i) {
        object_list[i] = it->second[0].as<std::string>();
      }
    }
    else {
      std::cout << "ERROR! Master Parameter file parameter " <<  it->first[0].as<std::string>() << " not recognized!\n";
      exit(1);
    }
  }
  write_objects(object_list);
  cleanup_files();
  delete[] object_list;
  return 0;
}

void init_files() {
  // Init files
  std::ofstream ppb_init("src/parse_params_body.cpp", std::ios_base::out);
  ppb_init << "// parse_params_body.cpp, generated automatically using make_params\n\nif (param_name.compare(\"n_runs\") == 0 || param_name.compare(\"run_name\") == 0) {}\n";
  ppb_init.close();
  std::ofstream prpb_init("src/print_params_body.cpp", std::ios_base::out);
  prpb_init << "// print_params_body.cpp, generated automatically using make_params\n\n";
  prpb_init.close();
  std::ofstream sp_init("src/parameters.h", std::ios_base::out);
  sp_init << "#ifndef _CYTOSCORE_PARAMETERS_H_\n#define _CYTOSCORE_PARAMETERS_H_\n\n// parameters.h, generated automatically using make_params\n\nstruct system_parameters {\n\n";
  sp_init.close();

}

void cleanup_files() {
  std::ofstream sp_cleanup("src/parameters.h", std::ios_base::out | std::ios_base::app);
  sp_cleanup << "\n};\n\n#endif // _CYTOSCORE_PARAMETERS_H_";
  sp_cleanup.close();
  std::ofstream ppb_cleanup("src/parse_params_body.cpp", std::ios_base::out | std::ios_base::app);
  ppb_cleanup << "else {\n  std::cout << \"  WARNING: Parameter \" << param_name << \" not recognized!\\n\";\n}\n";
  ppb_cleanup.close();
}

void write_params(std::string param_name, std::string param_value, std::string param_type) {
  write_parse_params_body_cpp(param_name, param_type);
  write_print_params_body_cpp(param_name);
  write_parameters_h(param_name, param_value, param_type);
}

void write_parameters_h(std::string param_name, std::string param_value, std::string param_type) {
  std::string file_name = "src/parameters.h";
  std::ofstream sp_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::string snippet = parameters_snippet(param_name, param_value, param_type);
  sp_file << snippet;
  sp_file.close();
}

void write_parse_params_body_cpp(std::string param_name, std::string param_type) {
  std::string file_name = "src/parse_params_body.cpp";
  std::ofstream ppb_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::string snippet = parse_params_snippet(param_name, param_type);
  ppb_file << snippet;
  ppb_file.close();
}

void write_print_params_body_cpp(std::string param_name) {
  std::string file_name = "src/print_params_body.cpp";
  std::ofstream prpb_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::ostringstream snippet;
  snippet << "param_file << \"" << param_name << " : \" << params." << param_name << " << \"\\n\";\n";
  prpb_file << snippet.str();
  prpb_file.close();
}

std::string parameters_snippet(std::string param_name, std::string param_value, std::string param_type) {

  std::ostringstream code;

  if (param_type.compare("string") == 0 || param_type.compare("char") == 0) {
    code << "  char *" << param_name << ";\n";
    return code.str();
  }
  else if (param_type.compare("int") == 0 ||
      param_type.compare("double") == 0 ||
      param_type.compare("long") == 0) {

    code << "  " << param_type << " " << param_name << " = " << param_value << ";\n";
    return code.str();
  }
  else {
    std::cout << "  ERROR! Parameter type not recognized.\n";
    exit(1);
  }
}

std::string parse_params_snippet(std::string param_name, std::string param_type) {

  std::ostringstream code;
  std::string atotype;
  if (param_type.compare("int") == 0)
    atotype = "atoi";
  else if (param_type.compare("double") == 0)
    atotype = "atof";
  else if (param_type.compare("long") == 0)
    atotype = "atol";
  else if (param_type.compare("string") == 0 || param_type.compare("char") == 0) {
    code <<  "else if ( param_name.compare(\"" << param_name << "\") == 0 ) {\n" <<
      "  params_[i_var]." << param_name << " = (char *) gmalloc((strlen(param_value.c_str()) + 1) * sizeof(char));\n"
      "  strcpy(params_[i_var]." << param_name << ", param_value.c_str());\n"
      "}\n";
    return code.str();
  }
  else  {
    std::cout << "  ERROR! Parameter type not recognized.\n";
    exit(1);
  }

  code << "else if (param_name.compare(\"" << param_name << "\") == 0){\n" <<
  "  params_[i_var]." << param_name << " = " << atotype << "(param_value.c_str());\n" <<
  //"  std::cout << \"  temperature = \" << param_value << std::endl;\n" <<
  "}\n";
  return code.str();
}

void write_objects(std::string *object_list) {

}

