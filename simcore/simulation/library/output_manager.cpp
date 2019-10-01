#include "output_manager.hpp"

void OutputManager::Init(system_parameters *params,
                         std::vector<SpeciesBase *> *species,
                         space_struct *space, int *i_step, std::string run_name,
                         bool reading_inputs, bool posits_only,
                         bool with_reloads, bool reduce_flag,
                         int reduce_factor) {
  run_name_ = run_name;
  reduce_flag_ = reduce_flag;
  reduce_factor_ = reduce_factor;
  with_reloads_ = with_reloads;
  species_ = species;
  params_ = params;
  space_ = space;
  i_step_ = i_step;
  n_posit_ = n_spec_ = n_checkpoint_ = ABS((int)params_->n_steps);
  posits_only_ = posits_only;
  n_thermo_ = params_->n_thermo;
  thermo_flag_ = params_->thermo_flag;
  if (!reading_inputs && thermo_flag_) {
    InitThermo(run_name_);
  } else if (reading_inputs && thermo_flag_ && thermo_analysis_) {
    InitThermoInput(run_name_);
  }
  std::string red_file_name =
      run_name_ + "_reduced" + std::to_string(reduce_factor);
  if (thermo_flag_ && reduce_flag_) {
    InitThermo(red_file_name);
  }

  for (auto it = species_->begin(); it != species_->end(); ++it) {
    if (!params_->load_checkpoint && reading_inputs) {
      (*it)->InitInputFiles(run_name_, posits_only_, with_reloads_);
      if (reduce_flag_) {
        (*it)->InitOutputFiles(red_file_name);
      } else if (params->checkpoint_from_spec) {
        (*it)->InitCheckpoints(run_name_);
      }
    } else {
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

void OutputManager::WriteOutputs() {
  if (posit_flag_ && (*i_step_ % n_posit_ == 0)) {
    WritePosits();
  }
  if (spec_flag_ && (*i_step_ % n_spec_ == 0)) {
    WriteSpecs();
  }
  if (checkpoint_flag_ && (*i_step_ % n_checkpoint_ == 0)) {
    WriteCheckpoints();
  }
  if (thermo_flag_ && (*i_step_ % n_thermo_ == 0)) {
    WriteThermo();
  }
}

void OutputManager::WritePosits() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() && *i_step_ % (*spec)->GetNPosit() == 0) {
      (*spec)->WritePosits();
    }
  }
}

void OutputManager::WriteSpecs() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetSpecFlag() && *i_step_ % (*spec)->GetNSpec() == 0) {
      (*spec)->WriteSpecs();
    }
  }
}

void OutputManager::WriteCheckpoints() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetCheckpointFlag() &&
        *i_step_ % (*spec)->GetNCheckpoint() == 0) {
      (*spec)->WriteCheckpoints();
    }
  }
}

void OutputManager::InitThermo(std::string fname) {
  fname.append(".thermo");
  othermo_file_.open(fname, std::ios::out | std::ios::binary);
  if (!othermo_file_.is_open()) {
    Logger::Error("Thermo file failed to open!");
  }
  othermo_file_.write(reinterpret_cast<char *>(&(params_->n_steps)),
                      sizeof(int));
  othermo_file_.write(reinterpret_cast<char *>(&(n_thermo_)), sizeof(int));
  othermo_file_.write(reinterpret_cast<char *>(&(params_->delta)),
                      sizeof(double));
  othermo_file_.write(reinterpret_cast<char *>(&(params_->n_dim)), sizeof(int));
}

void OutputManager::InitThermoInput(std::string fname) {
  fname.append(".thermo");
  ithermo_file_.open(fname, std::ios::in | std::ios::binary);
  if (!ithermo_file_.is_open()) {
    Logger::Error("Thermo file failed to open!");
  }
  int n_thermo, ndim, n_steps;
  double delta;
  ithermo_file_.read(reinterpret_cast<char *>(&n_steps), sizeof(int));
  ithermo_file_.read(reinterpret_cast<char *>(&n_thermo), sizeof(int));
  ithermo_file_.read(reinterpret_cast<char *>(&delta), sizeof(double));
  ithermo_file_.read(reinterpret_cast<char *>(&ndim), sizeof(int));
  if (n_steps != params_->n_steps || n_thermo != params_->n_thermo ||
      delta != params_->delta || ndim != params_->n_dim) {
    Logger::Error(
        "Input file %s does not match parameter file:\n"
        "n_steps: %d %d, n_thermo: %d %d, delta: %2.5f %2.5f, n_dim: %d %d",
        fname.c_str(), n_steps, params_->n_steps, n_thermo, params_->n_thermo, delta,
        params_->delta, ndim, params_->n_dim);
  }
}

void OutputManager::WriteThermo() {
  for (int i = 0; i < 9; ++i) {
    othermo_file_.write(reinterpret_cast<char *>(&(space_->unit_cell[i])),
                        sizeof(double));
  }
  for (int i = 0; i < 9; ++i) {
    othermo_file_.write(reinterpret_cast<char *>(&(space_->pressure_tensor[i])),
                        sizeof(double));
  }
  othermo_file_.write(reinterpret_cast<char *>(&(space_->pressure)),
                      sizeof(double));
  othermo_file_.write(reinterpret_cast<char *>(&(space_->volume)),
                      sizeof(double));
}

void OutputManager::ReadThermo() {
  for (int i = 0; i < 9; ++i) {
    ithermo_file_.read(reinterpret_cast<char *>(&(space_->unit_cell[i])),
                       sizeof(double));
  }
  for (int i = 0; i < 9; ++i) {
    ithermo_file_.read(reinterpret_cast<char *>(&(space_->pressure_tensor[i])),
                       sizeof(double));
  }
  ithermo_file_.read(reinterpret_cast<char *>(&(space_->pressure)),
                     sizeof(double));
  ithermo_file_.read(reinterpret_cast<char *>(&(space_->volume)),
                     sizeof(double));
}

void OutputManager::ReadInputs() {
  // if ( *i_step_ % n_posit_ != 0 && *i_step_ % n_spec_ != 0) {
  // return;
  //}
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if (posits_only_ && (*spec)->GetPositFlag() &&
        *i_step_ % (*spec)->GetNPosit() == 0) {
      (*spec)->ReadPosits();
    }
    // In the case that we only want to consider posits, but we only have spec
    // files, use specs to read average positions
    else if (posits_only_ && !(*spec)->GetPositFlag() &&
             (*spec)->GetSpecFlag() && *i_step_ % (*spec)->GetNSpec() == 0) {
      (*spec)->ReadPositsFromSpecs();
    } else if (!posits_only_ && (*spec)->GetSpecFlag() &&
               *i_step_ % (*spec)->GetNSpec() == 0) {
      (*spec)->ReadSpecs();
      if (params_->checkpoint_from_spec) {
        (*spec)->WriteCheckpoints();
      }
    }
  }
  if (thermo_flag_) {
    ReadThermo();
  }
  if (reduce_flag_) {
    WriteReduce();
  }
}

void OutputManager::WriteReduce() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() &&
        *i_step_ % (reduce_factor_ * (*spec)->GetNPosit()) == 0) {
      (*spec)->WritePosits();
    }
    if (!posits_only_ && (*spec)->GetSpecFlag() &&
        *i_step_ % (reduce_factor_ * (*spec)->GetNSpec()) == 0) {
      (*spec)->WriteSpecs();
    }
  }
  if (thermo_flag_ && *i_step_ % (reduce_factor_ * n_thermo_) == 0) {
    WriteThermo();
  }
}

void OutputManager::Close() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->CloseFiles();
  }
  if (thermo_flag_) {
    if (ithermo_file_.is_open()) {
      ithermo_file_.close();
    }
    if (othermo_file_.is_open()) {
      othermo_file_.close();
    }
  }
  if (params_->time_analysis) {
    WriteTime();
  }
}

void OutputManager::WriteTime() {
  std::string fname = run_name_;
  fname.append(".time");
  time_file_.open(fname, std::ios::out);
  time_file_ << "i_step n_steps delta\n";
  time_file_ << *i_step_ << " " << params_->n_steps << " " << params_->delta
             << "\n";
  time_file_.close();
}
