#include "species.h"

SpeciesBase::SpeciesBase(int n_members, system_parameters *params, space_struct *space, long seed) 
    : rng_(seed) {
  sid_ = SID::none;
  n_members_ = n_members;
  params_ = params;
  space_ = space;
  delta_ = params->delta;
}

void SpeciesBase::InitOutputFile(std::string run_name) {
  std::string sid_str = SIDToString(sid_);
  std::string file_name = run_name + "_" + sid_str + ".posit";
  std::cout<<"Posit file name is "<< run_name <<"_"<< sid_str <<" \n";
  oposit_file_.open(file_name, std::ios::out | std::ios::binary ); 
  if (!oposit_file_.is_open())
    std::cout<<"Output "<< file_name <<" file did not open\n";
  else{
    int size = sid_str.size();
    oposit_file_.write(reinterpret_cast<char*>(&size), sizeof(int));
    oposit_file_.write(sid_str.c_str(), sid_str.size());
    oposit_file_.write(reinterpret_cast<char*> (&params_->n_steps), sizeof(int));
    oposit_file_.write(reinterpret_cast<char*> (&params_->n_posit), sizeof(int));
  }
}

void SpeciesBase::InitInputFile(std::string in_file, std::ios::streampos beg){
  iposit_file_.open(in_file, std::ios::binary | std::ios::in );
  iposit_file_.seekg(beg);
}

void SpeciesBase::Close() { 
  if (oposit_file_.is_open())
    oposit_file_.close(); 
  if (iposit_file_.is_open())
    iposit_file_.close(); 
}

void SpeciesBase::InitConfig(system_parameters *params, space_struct *space, long seed) {
  n_members_ = 0;
  params_ = params;
  space_ = space;
  rng_.init(seed);
  delta_ = params->delta;
}

std::vector<Simple*> SpeciesBase::GetSimples() {
  std::vector<Simple*> sim;
  return sim;
}


