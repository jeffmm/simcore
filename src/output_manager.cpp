#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, int *i_step){
  params_ = params;
  i_step_ = i_step;
}

void OutputManager::WriteOutputs(){
  if ( (*i_step_ % params_->n_posit == 0) && 
       (params_->posit_flag == 1) ) 
    WriteSpeciesPosits();
}

void OutputManager::WriteSpeciesPosits(){
  for (auto spec_it : species_)
    spec_it.second->WritePosits(posit_files_[spec_it.first]);
}

void OutputManager::AddSpecie(SpeciesBase *spec){
  SID sid = spec->GetSID();
  species_[sid] = spec;
  if (params_->posit_flag == 1){
    std::string file_name = SIDToString(sid) + ".posit";
    posit_files_[sid];
    posit_files_[sid].open(file_name, std::ios::out | std::ios::binary);
    if (!posit_files_[sid].is_open())
      std::cout<<"Output "<< file_name <<" file did not open\n";
    }
}

std::ofstream& OutputManager::GetPositFile(SID sid){
  return posit_files_[sid];
}

