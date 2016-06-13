#include <iostream>
#include <yaml-cpp/yaml.h>
#include <string.h>
#include <sstream>
#include <ios>
#include <fstream>

void init_files();
void write_params(std::string param_name, std::string param_value, std::string param_type);
void write_parse_params_body_h(std::string param_name, std::string param_type);
void write_parameters_h(std::string param_name, std::string param_value, std::string param_type);
std::string parse_params_snippet(std::string param_name, std::string param_type);
std::string parameters_snippet(std::string param_name, std::string param_value, std::string param_type);
void write_print_params_body_h(std::string param_name);
void cleanup_files();
void write_objects(std::string *object_list,unsigned int size);

int main(int argc, char *argv[]) {

  std::string param_file, param_name, param_value, param_type;
  std::string *object_list;
  unsigned int size;
  if (argc==1) 
    param_file = "src/master_params.yaml";
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
      size = it->second.size();
      object_list = new std::string[size];
      for (unsigned int i=0; i<size; ++i) {
        object_list[i] = it->second[i].as<std::string>();
      }
    }
    else {
      std::cout << "ERROR! Master Parameter file parameter " <<  it->first[0].as<std::string>() << " not recognized!\n";
      exit(1);
    }
  }
  write_objects(object_list, size);
  cleanup_files();
  delete[] object_list;
  return 0;
}

void init_files() {
  // Init files
  std::ofstream ppb_init("src/parse_params_body.h", std::ios_base::out);
  ppb_init << "// parse_params_body.h, generated automatically using make_params\n\nif (param_name.compare(\"n_runs\") == 0 || param_name.compare(\"run_name\") == 0) {}\n";
  ppb_init.close();
  std::ofstream prpb_init("src/print_params_body.h", std::ios_base::out);
  prpb_init << "// print_params_body.h, generated automatically using make_params\n\n";
  prpb_init.close();
  std::ofstream sp_init("src/parameters.h", std::ios_base::out);
  sp_init << "#ifndef _SIMCORE_PARAMETERS_H_\n#define _SIMCORE_PARAMETERS_H_\n\n// parameters.h, generated automatically using make_params\n\nstruct system_parameters {\n\n";
  sp_init.close();

}

void cleanup_files() {
  std::ofstream sp_cleanup("src/parameters.h", std::ios_base::out | std::ios_base::app);
  sp_cleanup << "\n};\n\n#endif // _SIMCORE_PARAMETERS_H_";
  sp_cleanup.close();
  std::ofstream ppb_cleanup("src/parse_params_body.h", std::ios_base::out | std::ios_base::app);
  ppb_cleanup << "else {\n  std::cout << \"  WARNING: Parameter \" << param_name << \" not recognized!\\n\";\n}\n";
  ppb_cleanup.close();
}

void write_params(std::string param_name, std::string param_value, std::string param_type) {
  write_parse_params_body_h(param_name, param_type);
  write_print_params_body_h(param_name);
  write_parameters_h(param_name, param_value, param_type);
}

void write_parameters_h(std::string param_name, std::string param_value, std::string param_type) {
  std::string file_name = "src/parameters.h";
  std::ofstream sp_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::string snippet = parameters_snippet(param_name, param_value, param_type);
  sp_file << snippet;
  sp_file.close();
}

void write_parse_params_body_h(std::string param_name, std::string param_type) {
  std::string file_name = "src/parse_params_body.h";
  std::ofstream ppb_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::string snippet = parse_params_snippet(param_name, param_type);
  ppb_file << snippet;
  ppb_file.close();
}

void write_print_params_body_h(std::string param_name) {
  std::string file_name = "src/print_params_body.h";
  std::ofstream prpb_file(file_name.c_str(), std::ios_base::out | std::ios_base::app);
  std::ostringstream snippet;
  if (param_name.compare("grab_file")==0)
    snippet << "if (params.grab_flag) {\n  ";
  snippet << "param_file << \"" << param_name << " : \" << params." << param_name << " << \"\\n\";\n";
  if (param_name.compare("grab_file")==0)
    snippet << "}\n";
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

void write_objects(std::string *object_list, unsigned int size) {
  std::ostringstream code;
  code << "CC = g++-5\nCXX = g++-5\nCFLAGS = -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -g -DDEBUG -std=c++11\nLNKFLAGS = -lyaml-cpp -gnu\n\nifeq ($(CC),icpc)\n\tOMPFLAGS = -openmp -DBOB_OMP -Wno-deprecated\nelse\n\tCFLAGS += -Wno-deprecated-declarations -Wno-deprecated\n\tOMPFLAGS =	-fopenmp -DBOB_OMP\nendif\n\nGLXFLAGS = -L/opt/X11/lib -lglfw3 -framework OpenGL -lglew -lpthread\nGSLFLAGS = -L/usr/local/lib/ -I/usr/local/include/gsl/ -lgsl -lgslcblas\nFFTFLAGS = -L/usr/lib64 -lfftw3\nCXX = $(CC)\nCXXFLAGS = $(CFLAGS)\n\nSIMCORE_OBJS = simcore_main.o simulation_manager.o simulation.o object.o forces.o cell_list.o space.o allocate.o vector_algebra.o cpu.o error_exit.o generate_random_unit_vector.o graphics.o grabber.o writebmp.o minimum_distance.o\n\nCONFIGURE_OBJS = configure_simcore.o\n\nSIMCORE_OBJECT_OBJS =";
  for (unsigned int i=0; i<size; ++i)
    code << " " << object_list[i] << ".o";
  code << "\n\n../simcore:	$(SIMCORE_OBJS) $(SIMCORE_OBJECT_OBJS)\n\t$(CC)  -o ../simcore $(SIMCORE_OBJS) $(SIMCORE_OBJECT_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS) \n\n../configure_simcore: $(CONFIGURE_OBJS)\n\t$(CC) -o ../configure_simcore $(CONFIGURE_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS)\n";
  std::ofstream makefile("src/Makefile", std::ios_base::out);
  makefile << code.str();
  makefile.close();
  std::ostringstream header_code;
  for (unsigned int i=0; i<size; ++i)
    header_code << "#include \"" << object_list[i] << ".h\"\n";
  std::ofstream obj_header_file("src/objects.h", std::ios_base::out);
  obj_header_file << header_code.str();
  obj_header_file.close();
}

