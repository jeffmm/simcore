#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include <algorithm>

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cout << " Expected input file name. Exiting.\n";
    exit(1);
  }
  std::string posit_file_name = argv[1];
  int n_steps, n_posit;
  double delta, length, diameter, avg_pos[3], s_pos[3], avg_u[3];
  std::fstream iposit_file(posit_file_name, std::ios::in | std::ios::binary ); 
  std::fstream oposit_file(posit_file_name.append(".txt"), std::ios::out);
  if (!iposit_file.is_open()) {
    std::cout<<"ERROR: Input file "<< posit_file_name <<" did not open\n";
    exit(1);
  }
  iposit_file.read(reinterpret_cast<char*> (&n_steps), sizeof(int));
  iposit_file.read(reinterpret_cast<char*> (&n_posit), sizeof(int));
  iposit_file.read(reinterpret_cast<char*> (&delta), sizeof(double));
  oposit_file << "n_steps n_posit delta\n";
  oposit_file << n_steps << " " << n_posit << " " << delta << "\n";
  int n_points = (int) floor(n_steps / n_posit);
  oposit_file << "position(3) \n"; //scaled_position(3) orientation(3) diameter length\n";
  for (int i=0; i<n_points; ++i) {
    if (iposit_file.eof()) break;
    int size;
    iposit_file.read(reinterpret_cast<char*> (&size), sizeof(size));
    double avg_pos[3], avg_u[3], s_pos[3];
    for (int i=0; i<3; ++i)
      iposit_file.read(reinterpret_cast<char*>(&avg_pos[i]), sizeof(double));
    for (int i=0; i<3; ++i)
      iposit_file.read(reinterpret_cast<char*>(&s_pos[i]), sizeof(double));
    for (int i=0; i<3; ++i)
      iposit_file.read(reinterpret_cast<char*>(&avg_u[i]), sizeof(double));
    iposit_file.read(reinterpret_cast<char*>(&diameter), sizeof(diameter));
    iposit_file.read(reinterpret_cast<char*>(&length), sizeof(length));
    oposit_file << avg_pos[0]<< " " <<avg_pos[1] << "\n"; // << " " <<avg_pos[2]<< " " <<
        //s_pos[0]<< " " << s_pos[1]<< " " << s_pos[2]<< " " <<
        //avg_u[0]<< " " <<avg_u[1]<< " " <<avg_u[2]<< " " <<
        //diameter<< " " << length << "\n";
  }
  iposit_file.close();
  oposit_file.close();
  return 0;
}
