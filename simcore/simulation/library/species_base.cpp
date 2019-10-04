#include "species_base.hpp"

void SpeciesBase::Init(system_parameters *params, species_base_parameters *sparams, space_struct *space) {
  params_ = params;
  space_ = space;
  n_members_ = 0;
  Logger::Debug("Initializing species %s", GetSID()._to_string());
}

void SpeciesBase::InitPositFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name = run_name + "_" + sid_str + "_" + spec_label_ + ".posit";
  oposit_file_.open(posit_file_name, std::ios::out | std::ios::binary);
  if (!oposit_file_.is_open()) {
    Logger::Error("Output file %s did not open", posit_file_name.c_str());
  }
  int n_posit = GetNPosit();
  oposit_file_.write(reinterpret_cast<char *>(&params_->n_steps), sizeof(int));
  oposit_file_.write(reinterpret_cast<char *>(&n_posit), sizeof(int));
  oposit_file_.write(reinterpret_cast<char *>(&params_->delta), sizeof(double));
}

void SpeciesBase::InitPositFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name = run_name + "_" + sid_str + "_" + spec_label_ + ".posit";
  iposit_file_.open(posit_file_name, std::ios::in | std::ios::binary);
  if (!iposit_file_.is_open()) {
    Logger::Error("Input file %s did not open", posit_file_name.c_str());
  }
  // long n_steps;
  int n_posit, n_steps;
  double delta;
  iposit_file_.read(reinterpret_cast<char *>(&n_steps), sizeof(int));
  iposit_file_.read(reinterpret_cast<char *>(&n_posit), sizeof(int));
  iposit_file_.read(reinterpret_cast<char *>(&delta), sizeof(double));
  if (n_steps != params_->n_steps || n_posit != GetNPosit() ||
      delta != params_->delta) {
    Logger::Warning("Input file %s does not match parameter file\n",
            posit_file_name.c_str());
    printf("%d %d, %d %d, %2.5f %2.5f\n", n_steps, params_->n_steps, n_posit,
           GetNPosit(), delta, params_->delta);
  }
  ReadPosits();
}

void SpeciesBase::InitSpecFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name = run_name + "_" + sid_str + "_" + spec_label_ + ".spec";
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary);
  if (!ospec_file_.is_open()) {
    Logger::Error("Output file %s did not open", spec_file_name.c_str());
  }
  int n_spec = GetNSpec();
  ospec_file_.write(reinterpret_cast<char *>(&params_->n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&n_spec), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&params_->delta), sizeof(double));
}

bool SpeciesBase::HandleEOF() {
  if (++spec_file_iterator_) {
    if (params_->i_step < params_->n_steps_equil) {
      params_->n_steps_equil -= params_->i_step;
    } else {
      params_->n_steps_equil = 0;
    }
    params_->i_step = 0;
    std::ostringstream file_name, nload;
    std::string sid_str = sid_._to_string();
    file_name << params_->run_name;
    nload << std::setw(3) << std::setfill('0') << spec_file_iterator_;
    size_t pos;
    if ((pos = file_name.str().find("reload")) == std::string::npos) {
      // This is not a reload file currently
      if (params_->reduced) {
        /* The file is either a reduced file, or we are currently reducing */

        if ((pos = file_name.str().find("_reduced")) == std::string::npos) {
          /* we are currently reducing, so input file does not have reduce in
           * name */
          pos = file_name.tellp();
          file_name << "_reload" << nload.str();
        } else {
          if (!params_->reload_reduce_switch) {
            /* need to (probably) prefix with reload, assuming the reduction
               came after the reload (most typical case) */
            file_name.seekp(pos);
            file_name << "_reload" << nload.str();
            file_name << "_reduced" << params_->reduced;
          } else {
            file_name << "_reload" << nload.str();
          }
        }
      } else {
        file_name << "_reload" << nload.str();
      }
    } else {
      // The file is currently a reload file, simply seek to beginning of
      // substring "reload"
      file_name.seekp(pos);
      file_name << "_reload" << nload.str();
      if (params_->reduced) {
        file_name << "_reduced" << params_->reduced;
      }
    }
    file_name << "_" << sid_str << "_" << spec_label_ << ".spec";
    if (ispec_file_.is_open()) {
      ispec_file_.close();
    }
    printf("\n  Initializing new input file, %s\n", file_name.str().c_str());
    return InitSpecFileInputFromFile(file_name.str());
  } else {
    return false;
  }
}

bool SpeciesBase::InitSpecFileInputFromFile(std::string spec_file_name) {
  ispec_file_.open(spec_file_name, std::ios::in | std::ios::binary);
  if (!ispec_file_.is_open()) {
    return false;
  }
  // long n_steps;
  int n_spec, n_steps;
  double delta;
  ispec_file_.read(reinterpret_cast<char *>(&n_steps), sizeof(int));
  ispec_file_.read(reinterpret_cast<char *>(&n_spec), sizeof(int));
  ispec_file_.read(reinterpret_cast<char *>(&delta), sizeof(double));
  if (n_steps != params_->n_steps || n_spec != GetNSpec() ||
      delta != params_->delta) {
    printf("\nn_steps: %d %d, ", n_steps, params_->n_steps);
    printf("n_spec: %d %d, ", n_spec, GetNSpec());
    printf("delta: %2.2f %2.2f\n", delta, params_->delta);
    Logger::Warning("Input file %s does not match parameter file\n",
            spec_file_name.c_str());
  }
  ReadSpecs();
  return true;
}

void SpeciesBase::InitSpecFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name = run_name + "_" + sid_str + "_" + spec_label_ + ".spec";
  if (!InitSpecFileInputFromFile(spec_file_name)) {
    Logger::Error("Input file %s did not open", spec_file_name.c_str());
  }
}

void SpeciesBase::InitOutputFiles(std::string run_name) {
  if (GetPositFlag())
    InitPositFile(run_name);
  if (GetSpecFlag())
    InitSpecFile(run_name);
  if (GetCheckpointFlag())
    InitCheckpoints(run_name);
}

void SpeciesBase::InitCheckpoints(std::string run_name) {
  std::string sid_str = sid_._to_string();
  checkpoint_file_ = run_name + "_" + sid_str + "_" + spec_label_ + ".checkpoint";
}

void SpeciesBase::LoadFromCheckpoints(std::string run_name,
                                      std::string checkpoint_run_name) {
  std::string sid_str = sid_._to_string();
  checkpoint_file_ = checkpoint_run_name + "_" + sid_str + "_" + spec_label_ + ".checkpoint";
  if (!GetCheckpointFlag()) {
    Logger::Error("Checkpoint file %s not available for parameter file!",
                  checkpoint_file_.c_str());
  }
  ReadCheckpoints();
  InitOutputFiles(run_name);
}

void SpeciesBase::InitInputFiles(std::string run_name, bool posits_only,
                                 bool with_reloads) {
  if (posits_only && GetPositFlag()) {
    InitPositFileInput(run_name);
  } else if (GetSpecFlag() && with_reloads) {
    spec_file_iterator_ = 0;
    InitSpecFileInput(run_name);
  } else if (GetSpecFlag()) {
    InitSpecFileInput(run_name);
  }
}

void SpeciesBase::CloseFiles() {
  if (oposit_file_.is_open())
    oposit_file_.close();
  if (iposit_file_.is_open())
    iposit_file_.close();
  if (ospec_file_.is_open())
    ospec_file_.close();
  if (ispec_file_.is_open())
    ispec_file_.close();
  // FinalizeAnalysis();
}
