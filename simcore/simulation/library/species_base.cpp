#include "species_base.hpp"


SpeciesBase::SpeciesBase(unsigned long seed) : rng_(seed) {}
void SpeciesBase::SetSID(species_id sid) { sid_ = sid; }

/* Static functions*/
void SpeciesBase::SetParams(system_parameters *params) { params_ = params; }
void SpeciesBase::SetSpace(space_struct *space) { space_ = space; }
const space_struct *SpeciesBase::space_ = nullptr;
const system_parameters *SpeciesBase::params_ = nullptr;;

void SpeciesBase::InitPositFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name =
      run_name + "_" + sid_str + "_" + GetSpeciesName() + ".posit";
  oposit_file_.open(posit_file_name, std::ios::out | std::ios::binary);
  if (!oposit_file_.is_open()) {
    Logger::Error("Output file %s did not open", posit_file_name.c_str());
  }
  int n_posit = GetNPosit();
  int n_steps = params_->n_steps;
  double delta = params_->delta;
  oposit_file_.write(reinterpret_cast<char *>(&n_steps), sizeof(int));
  oposit_file_.write(reinterpret_cast<char *>(&n_posit), sizeof(int));
  oposit_file_.write(reinterpret_cast<char *>(&delta), sizeof(double));
}

void SpeciesBase::InitPositFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string posit_file_name =
      run_name + "_" + sid_str + "_" + GetSpeciesName() + ".posit";
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
    Logger::Warning("Input file %s does not match parameter file: "
                    "n_steps: %d %d, "
                    "n_spec: %d %d, "
                    "delta: %2.2f %2.2f ",
                    posit_file_name.c_str(), n_steps, params_->n_steps, n_posit,
                    GetNPosit(), delta, params_->delta);
  }
  ReadPosits();
}

void SpeciesBase::InitSpecFile(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name =
      run_name + "_" + sid_str + "_" + GetSpeciesName() + ".spec";
  Logger::Trace("Initializing spec file %s", spec_file_name.c_str());
  ospec_file_.open(spec_file_name, std::ios::out | std::ios::binary);
  if (!ospec_file_.is_open()) {
    Logger::Error("Output file %s did not open", spec_file_name.c_str());
  }
  int n_spec = GetNSpec();
  int n_steps = params_->n_steps;
  double delta = params_->delta;
  ospec_file_.write(reinterpret_cast<char *>(&n_steps), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&n_spec), sizeof(int));
  ospec_file_.write(reinterpret_cast<char *>(&delta), sizeof(double));
}

bool SpeciesBase::HandleEOF() {
  if (++spec_file_iterator_) {
    //if (params_->i_step < params_->n_steps_equil) {
      //params_->n_steps_equil -= params_->i_step;
    //} else {
      //params_->n_steps_equil = 0;
    //}
    //params_->i_step = 0;
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
    file_name << "_" << sid_str << "_" << GetSpeciesName() << ".spec";
    if (ispec_file_.is_open()) {
      ispec_file_.close();
    }
    Logger::Info("Switching to new spec file, %s", file_name.str().c_str());
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
    Logger::Warning("Input file %s does not match parameter file: "
                    "n_steps: %d %d, "
                    "n_spec: %d %d, "
                    "delta: %2.2f %2.2f ",
                    spec_file_name.c_str(), n_steps, params_->n_steps, n_spec,
                    GetNSpec(), delta, params_->delta);
  }
  ReadSpecs();
  return true;
}

void SpeciesBase::InitSpecFileInput(std::string run_name) {
  std::string sid_str = sid_._to_string();
  std::string spec_file_name =
      run_name + "_" + sid_str + "_" + GetSpeciesName() + ".spec";
  if (!InitSpecFileInputFromFile(spec_file_name)) {
    Logger::Error("Input file %s did not open", spec_file_name.c_str());
  }
}

void SpeciesBase::InitOutputFiles(std::string run_name) {
  Logger::Trace("Initializing output files for %s %s", sid_._to_string(),
                GetSpeciesName().c_str());
  if (GetPositFlag())
    InitPositFile(run_name);
  if (GetSpecFlag())
    InitSpecFile(run_name);
  if (GetCheckpointFlag())
    InitCheckpoints(run_name);
}

void SpeciesBase::InitCheckpoints(std::string run_name) {
  std::string sid_str = sid_._to_string();
  checkpoint_file_ =
      run_name + "_" + sid_str + "_" + GetSpeciesName() + ".checkpoint";
}

void SpeciesBase::LoadFromCheckpoints(std::string run_name,
                                      std::string checkpoint_run_name) {
  std::string sid_str = sid_._to_string();
  checkpoint_file_ = checkpoint_run_name + "_" + sid_str + "_" +
                     GetSpeciesName() + ".checkpoint";
  Logger::Trace("Loading %s %s from checkpoint file %s", sid_._to_string(),
                GetSpeciesName().c_str(), checkpoint_file_.c_str());
  if (!GetCheckpointFlag()) {
    Logger::Error("Checkpoint file %s not available for parameter file!",
                  checkpoint_file_.c_str());
  }
  ReadCheckpoints();
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
  Logger::Trace("Closing output files for %s %s", sid_._to_string(),
                GetSpeciesName().c_str());
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
