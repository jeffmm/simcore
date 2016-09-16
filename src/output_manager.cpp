#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, int *i_step){
  params_ = params;
  i_step_ = i_step;
}

void OutputManager::WriteOutputs(){

  if ( (*i_step_ % params_->n_posit == 0) && params_->posit_flag)  {
    WriteSpeciesPosits();
  }
}

void OutputManager::WriteSpeciesPosits(){
  for (auto spec_it : species_)
    spec_it.second->WritePosits();
}

void OutputManager::AddSpecie(SpeciesBase *spec){
  SID sid = spec->GetSID();
  species_[sid] = spec;
  if (params_->posit_flag == 1){
    spec->InitOutputFile();
  }
}

//std::ofstream& OutputManager::GetPositFile(SID sid){
  //return oposit_files_[sid];
//}

void OutputManager::Close() {
  for (auto spec : species_){
    spec.second->Close();
  }
}

void OutputManager::Clear() {}

