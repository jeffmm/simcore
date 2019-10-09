#ifndef _SIMCORE_OUTPUT_MANAGER_H_
#define _SIMCORE_OUTPUT_MANAGER_H_

#include "auxiliary.hpp"
#include "object.hpp"
#include "species.hpp"
#include "parse_flags.hpp"

class OutputManager {
private:
  bool posit_flag_ = false;
  bool spec_flag_ = false;
  bool checkpoint_flag_ = false;
  bool thermo_flag_ = false;
  bool posits_only_ = false;
  bool reduce_flag_ = false;
  bool thermo_analysis_ = false;
  bool with_reloads_ = false;
  int n_posit_;
  int n_spec_;
  int n_checkpoint_;
  int n_thermo_;
  int reduce_factor_;
  std::string run_name_;
  system_parameters *params_;
  std::vector<SpeciesBase *> *species_;
  space_struct *space_;
  void WritePosits();
  void WriteSpecs();
  void WriteCheckpoints();
  void WriteThermo();
  void ReadThermo();
  void InitThermo(std::string fname);
  void InitThermoInput(std::string fname);
  void WriteReduce();
  void WriteTime();
  std::fstream othermo_file_;
  std::fstream ithermo_file_;
  std::fstream time_file_;

public:
  OutputManager() {}
  void Init(system_parameters *params, std::vector<SpeciesBase *> *species,
            space_struct *space, bool reading_inputs = false,
            run_options *run_opts = nullptr);
  // void SetPosits(std::vector<std::string> posit_files) {
  // posit_files_ = posit_files;
  //}
  int GetNPosit() { return n_posit_; }
  int GetNSpec() { return n_spec_; }
  int GetNCheckpoint() { return n_checkpoint_; }
  void WriteOutputs();
  void InitInputs();
  void ReadInputs();
  void Close();
};

#endif //_SIMCORE_OUTPUT_MANAGER_H_
