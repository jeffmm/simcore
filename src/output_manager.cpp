#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, 
                         std::vector<SpeciesBase*> *species,
                         int *i_step, std::string run_name) {
  run_name_ = run_name;
  n_posit_ = params->n_posit;
  species_ = species;
  params_ = params;
  i_step_ = i_step;
}

void OutputManager::WriteOutputs(){
  if ( (*i_step_ % n_posit_ == 0) && params_->posit_flag) 
    WritePosits();
}

void OutputManager::WritePosits(){
  //for (auto spec_it : species_)
    //spec_it->WritePosits();
}

void OutputManager::Close() {
  //for (auto spec : species_){
    //spec.second->Close();
  //}
}

void OutputManager::InitPositInput() {
  //for (auto pos_it : posit_files_){
    //int nchar, n_steps;
    //std::fstream ip;

    ////Open binary file of all the positions
    //ip.open(pos_it, std::ios::binary | std::ios::in );
    //if (!ip.is_open())
      //error_exit("ERROR: Input file %s did not open.\n",pos_it.c_str());
    ////Get Sim data from posit file
    //ip.read(reinterpret_cast<char*>(&nchar), sizeof(int));
    //std::string sid_str(nchar, ' ');
    //ip.read(&sid_str[0], nchar);

    //SID sid_posit = StringToSID(sid_str);

    //ip.read(reinterpret_cast<char*>(&(n_steps)), sizeof(int));
    //ip.read(reinterpret_cast<char*>(&(n_posit_)), sizeof(int));

    ////Need n_steps to be the lowest
    //if (params_->n_steps > n_steps) params_->n_steps = n_steps;

    ////Get the location where species data starts and close input for now
    //std::ios::streampos beg = ip.tellg();
    //ip.close();

    ////Find species that correlates to SID and add it to specs
    //for (auto spec_it : species_ )
      //if (sid_posit == spec_it.first){
        //spec_it.second->InitInputFile(pos_it, beg);
        //specs[sid_posit] = spec_it.second;
        //break;
      //}
  //}
}

void OutputManager::ReadPosits() {
  //for( auto ps_it : posit_spec_) 
    //if (*i_step_ % ps_it.second == 0 )
      //species_[ps_it.first]->ReadPosits();
}

