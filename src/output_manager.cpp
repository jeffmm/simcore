#include "output_manager.h"

OutputManager::OutputManager(){}

void OutputManager::Init(system_parameters *params, 
                         std::vector<SpeciesBase*> *species,
                         int *i_step, std::string run_name) {
  run_name_ = run_name;
  species_ = species;
  params_ = params;
  i_step_ = i_step;
  posit_flag_ = false;
  // Set maximum possible n_posit, then determine smallest species n_posit
  n_posit_ = params->n_steps; 
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag()) {
      posit_flag_ = true;
      (*spec)->InitOutputFile(run_name_);
    }
    if ((*spec)->GetNPosit() < n_posit_)
      n_posit_ = (*spec)->GetNPosit();
  }
  if (make_movie_) {
    InitPositInput();
    ReadPosits();
  }
}

void OutputManager::WriteOutputs(){
  if ( posit_flag_ && (*i_step_ % n_posit_ == 0) ) 
    WritePosits();
}

void OutputManager::WritePosits(){
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ( (*spec)->GetPositFlag() && *i_step_ % (*spec)->GetNPosit() == 0 ) 
      (*spec)->WritePosits();
  }
}

void OutputManager::ReadPosits() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() && *i_step_ % (*spec)->GetNPosit() == 0 )
      (*spec)->ReadPosits();
  }
}

void OutputManager::Close() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->ClosePosit();
  }
}

void OutputManager::InitPositInput() {
  for (auto pos_it=posit_files_.begin(); pos_it!=posit_files_.end(); ++pos_it){
    int nchar, n_steps, n_pos;
    std::fstream ip;

    //Open binary file of all the positions
    ip.open(*pos_it, std::ios::binary | std::ios::in );
    if (!ip.is_open())
      error_exit("ERROR: Input file %s did not open.\n",pos_it->c_str());
    //Get Sim data from posit file
    ip.read(reinterpret_cast<char*>(&nchar), sizeof(int));
    std::string sid_str(nchar, ' ');
    ip.read(&sid_str[0], nchar);

    SID sid_posit = StringToSID(sid_str);

    ip.read(reinterpret_cast<char*>(&(n_steps)), sizeof(int));
    ip.read(reinterpret_cast<char*>(&(n_pos)), sizeof(int));

    //Need n_steps and n_posit_ to be the lowest
    if (params_->n_steps > n_steps)
      params_->n_steps = n_steps;
    if (n_posit_ > n_pos) 
      n_posit_ = n_pos;

    //Get the location where species data starts and close input for now
    std::ios::streampos beg = ip.tellg();
    ip.close();

    //Find species that correlates to SID and add it to specs
    for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
      if (sid_posit == (*spec)->GetSID()){
        (*spec)->SetPositFlag(true);
        (*spec)->SetNPosit(n_pos);
        (*spec)->InitInputFile(*pos_it, beg);
        break;
      }
    }
  }
}

