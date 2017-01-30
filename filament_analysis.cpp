#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <signal.h>
#include <string.h>
#include "math.h"
#include <ios>
#include <fstream>
#include <sstream>
#include <yaml-cpp/yaml.h>
  
bool debug = true; // debug mode
void print(const char *print_msg, ...);

int main(int argc, char *argv[]) {

  // Catch usage error
  if (argc != 3) {
    std::cout << "Usage: \n  " << argv[0] << " param_file posit_file\n";
    exit(1);
  }

  std::fstream ip;
  std::string param_file = argv[1];
  std::string posit_file = argv[2];
  YAML::Node params = YAML::LoadFile(param_file);
  double delta = params["delta"].as<double>();
  double persistence_length = 80;
  //double persistence_length = params["persistence_length"].as<double>();
  // open posit file
  {
    ip.open(posit_file.c_str(), std::ios::binary | std::ios::in);
    if (!ip.is_open()) {
      std::cout << "Failed to open posit file " << posit_file << "!\n";
      exit(1);
    }
  }

  // read header file
  int n_char, n_steps, n_posit;
  {
    ip.read(reinterpret_cast<char*>(&n_char), sizeof(int));
    std::string sid_str(n_char, ' ');
    ip.read(&sid_str[0], n_char);
    ip.read(reinterpret_cast<char*>(&(n_steps)), sizeof(int));
    ip.read(reinterpret_cast<char*>(&(n_posit)), sizeof(int));
    print(sid_str.c_str());
    print("\n");
    print("n_steps: %d\n", n_steps);
    print("n_posit: %d\n", n_posit);
  }

  int n_iterations = n_steps / n_posit;
  int n_objs, n_objs_alloc, n_bonds, n_bonds_alloc;
  double diameter, length, child_length,
         **head_pos, **tail_pos, ***u_bond;
  {
    ip.read(reinterpret_cast<char*>(&n_objs), sizeof(n_objs));
    for (int i_obj=0; i_obj<n_objs; ++i_obj) {
      double temp,temp_array[3];
      for (int i=0; i<3; ++i)
        ip.read(reinterpret_cast<char*>(&temp), sizeof(double));
      ip.read(reinterpret_cast<char*>(&n_bonds), sizeof(int));
      for (int i=0; i<n_bonds+2; ++i) 
        ip.read(reinterpret_cast<char*>(&temp_array), 3*sizeof(double));
    }
    // Allocate arrays
    print("n_objs: %d\n",n_objs);
    print("n_bonds: %d\n", n_bonds);
    n_objs_alloc = n_objs;
    n_bonds_alloc = n_bonds;
    head_pos = new double*[n_objs];
    tail_pos = new double*[n_objs];
    u_bond = new double**[n_objs];
    for (int i=0; i<n_objs; ++i) {
      head_pos[i] = new double[3];
      tail_pos[i] = new double[3];
      u_bond[i] = new double*[n_bonds];
      for (int j=0; j<n_bonds; ++j) {
        u_bond[i][j] = new double[3];
      }
    }
  }

  // Loop through posit data
  double avg_e2e_sqr, avg_e2e, avg_bend_energy;
  avg_e2e_sqr = avg_e2e = avg_bend_energy = 0;
  {
    for (int i_iteration=0; i_iteration<n_iterations; ++i_iteration) {
      if (ip.eof()) 
        break;
      ip.read(reinterpret_cast<char*>(&n_objs), sizeof(n_objs));
      if (n_objs != n_objs_alloc) {
        printf("Error! Filament analysis does not support differing number of objects.\n");
        exit(1);
      }
      for (int i_obj=0; i_obj<n_objs; ++i_obj) {
        ip.read(reinterpret_cast<char*>(&diameter), sizeof(double));
        ip.read(reinterpret_cast<char*>(&length), sizeof(double));
        ip.read(reinterpret_cast<char*>(&child_length), sizeof(double));
        ip.read(reinterpret_cast<char*>(&n_bonds), sizeof(int));
        if (n_bonds != n_bonds_alloc) {
          printf("Error! Filament analysis does not support different number of bonds.\n");
          exit(1);
        }
        for (int i=0; i<3; ++i)
          ip.read(reinterpret_cast<char*>(&head_pos[i_obj][i]), sizeof(double));
        for (int i=0; i<3; ++i)
          ip.read(reinterpret_cast<char*>(&tail_pos[i_obj][i]), sizeof(double));
        for (int i_bond=0; i_bond<n_bonds; ++i_bond) {
          for (int i=0; i<3; ++i)
            ip.read(reinterpret_cast<char*>(&u_bond[i_obj][i_bond][i]), sizeof(double));
        }
      }


      {
        double cos_angle, bend_energy, e2edist;
        double bend_energy_zero = -(n_bonds-1)*persistence_length/child_length;
        for (int i_obj=0; i_obj<n_objs; ++i_obj) {
          e2edist = bend_energy = 0;
          for (int i=0; i<3; ++i) {
            double diff = tail_pos[i_obj][i] - head_pos[i_obj][i];
            e2edist += diff*diff;
          }
          avg_e2e_sqr += e2edist;
          e2edist = sqrt(e2edist);
          avg_e2e += e2edist;
          for (int i_bond=0; i_bond<n_bonds-1; ++i_bond) {
            cos_angle = 0;
            for (int i=0; i<3; ++i)
              cos_angle += u_bond[i_obj][i_bond][i]*u_bond[i_obj][i_bond+1][i];
            //printf("cos_angle: %2.8f\n", cos_angle);
            bend_energy += cos_angle;
          }
          bend_energy*= -persistence_length/child_length;
          bend_energy -= bend_energy_zero;
          avg_bend_energy += bend_energy;
          //print("obj %d:\n",i_obj);
          //print("  bend_energy: %2.8f\n", bend_energy);
          //print("  e2edist: %2.8f\n", e2edist);
          //exit(0);
        }

      }
    }
  }

  avg_bend_energy /= n_iterations;
  avg_e2e_sqr /= n_iterations;
  avg_e2e /= n_iterations;
  printf("avg bend_energy: %2.8f\n",avg_bend_energy);
  printf("avg end-to-end dist: %2.8f\n",avg_e2e);
  printf("avg sqr end-to-end dist: %2.8f\n",avg_e2e_sqr);
  printf("fluctuations: %2.8f\n",sqrt(avg_e2e_sqr-avg_e2e*avg_e2e));

  ip.close();
  exit(0);
}

void print(const char *print_msg, ...) {
  va_list args;
  va_start(args, print_msg);
  vfprintf(stdout, print_msg, args);
  va_end(args);
  return;
}

