#include "output_manager.h"

void OutputManager::Init(system_parameters *params, 
                         std::vector<SpeciesBase*> *species,
                         space_struct *space,
                         int *i_step, std::string run_name,
                         bool reading_inputs, bool posits_only) {
  run_name_ = run_name;
  species_ = species;
  params_ = params;
  space_ = space;
  i_step_ = i_step;
  n_posit_ = n_spec_ = n_checkpoint_ = params_->n_steps;
  posits_only_ = posits_only;
  if (thermo_flag_ = params_->thermo_flag) {
    n_thermo_ = params_->n_thermo;
    InitThermo();
  }
  for (auto it = species_->begin(); it != species_->end(); ++it) {
    if (params->load_checkpoint) {
      (*it)->InitCheckpoints(run_name_);
    }
    else if (reading_inputs) {
      (*it)->InitInputFiles(run_name_, posits_only);
    }
    else {
      (*it)->InitOutputFiles(run_name_);
    }
    if ((*it)->GetPositFlag()) {
      posit_flag_ = true;
      if ((*it)->GetNPosit() < n_posit_) {
        n_posit_ = (*it)->GetNPosit();
      }
    }
    if ((*it)->GetSpecFlag()) {
      spec_flag_ = true;
      if ((*it)->GetNSpec() < n_spec_) {
        n_spec_ = (*it)->GetNSpec();
      }
    }
    if ((*it)->GetCheckpointFlag()) {
      checkpoint_flag_ = true;
      if ((*it)->GetNCheckpoint() < n_checkpoint_) {
        n_checkpoint_ = (*it)->GetNCheckpoint();
      }
    }
  }
}

void OutputManager::WriteOutputs(){
  if ( posit_flag_ && (*i_step_ % n_posit_ == 0) ) {
    WritePosits();
  }
  if ( spec_flag_ && (*i_step_ % n_spec_ == 0) ) {
    WriteSpecs();
  }
  if ( checkpoint_flag_ && (*i_step_ % n_checkpoint_ == 0) ) {
    WriteCheckpoints();
  }
  if ( thermo_flag_ && (*i_step_ % n_thermo_ == 0) ) {
    WriteThermo();
  }
}

void OutputManager::WritePosits() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ( (*spec)->GetPositFlag() && *i_step_ % (*spec)->GetNPosit() == 0 ) {
      (*spec)->WritePosits();
    }
  }
}

void OutputManager::WriteSpecs() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ( (*spec)->GetSpecFlag() && *i_step_ % (*spec)->GetNSpec() == 0 ) {
      (*spec)->WriteSpecs();
    }
  }
}

void OutputManager::WriteCheckpoints() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ( (*spec)->GetCheckpointFlag() && *i_step_ % (*spec)->GetNCheckpoint() == 0 ) {
      (*spec)->WriteCheckpoints();
    }
  }
}

void OutputManager::InitThermo() {
  std::string fname = run_name_;
  fname.append(".thermo");
  thermo_file_.open(fname, std::ios::out | std::ios::binary);
  if (!thermo_file_.is_open()) {
    std::cout << "ERROR: Thermo file failed to open!\n";
    exit(1);
  }
  thermo_file_.write(reinterpret_cast<char*>(&(params_->n_steps)), sizeof(int));
  thermo_file_.write(reinterpret_cast<char*>(&(n_thermo_)), sizeof(int));
  thermo_file_.write(reinterpret_cast<char*>(&(params_->delta)), sizeof(double));
  thermo_file_.write(reinterpret_cast<char*>(&(params_->n_dim)), sizeof(int));
}

void OutputManager::WriteThermo() {
  for (int i=0; i<9; ++i) {
    thermo_file_.write(reinterpret_cast<char*>(&(space_->unit_cell[i])), sizeof(double));
  }
  for (int i=0; i<9; ++i) {
    thermo_file_.write(reinterpret_cast<char*>(&(space_->pressure_tensor[i])), sizeof(double));
  }
  thermo_file_.write(reinterpret_cast<char*>(&(space_->pressure)), sizeof(double));
  thermo_file_.write(reinterpret_cast<char*>(&(space_->volume)), sizeof(double));
}

void OutputManager::ReadInputs() {
  if ( *i_step_ % n_posit_ != 0 && *i_step_ % n_spec_ != 0) {
    return;
  }
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() && *i_step_ % (*spec)->GetNPosit() == 0 ) {
      (*spec)->ReadPosits();
    }
    if (!posits_only_ && (*spec)->GetSpecFlag() && *i_step_ % (*spec)->GetNSpec() == 0 ) {
      (*spec)->ReadSpecs();
    }
  }
}

void OutputManager::Close() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->CloseFiles();
  }
  if (thermo_flag_ && thermo_file_.is_open()) {
    thermo_file_.close();
  }
}

