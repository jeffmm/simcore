#ifndef _SIMCORE_ANALYSIS_H_
#define _SIMCORE_ANALYSIS_H_

#include "auxiliary.hpp"

class AnalysisBase {
private:
  std::string analysis_name_ = "unnamed";
  bool run_interaction_analysis_ = false;

protected:
  std::fstream output_;
  static const system_parameters *params_;
  static const space_struct *space_;
  void SetAnalysisName(std::string name) { analysis_name_ = name; }
  std::string GetAnalysisName() { return analysis_name_; }
  void RequireInteractionAnalysis() { run_interaction_analysis_ = true; }

public:
  static void SetParams(system_parameters *params);
  static void SetSpace(space_struct *space);
  virtual void Init() {
    if (params_ == nullptr || space_ == nullptr) {
      Logger::Error("System parameters or space data structure are unset in "
                    "AnalysisBase");
    }
  }
  virtual void Run() {}
  virtual void End() {
    if (output_.is_open()) {
      output_.close();
    }
  }
  bool CheckInteractionAnalysis() { return run_interaction_analysis_; }
  virtual void Draw(std::vector<graph_struct *> &graph_array) {}
};

template <typename T, unsigned char S> class Analysis : public AnalysisBase {
protected:
  std::vector<T> *members_;
  species_parameters<S> *sparams_;
  species_id sid_;
  int n_members_ = -1;
  double time_ = 0;
  int iteration_ = 0;

  virtual void InitOutput() {
    std::string file_name = params_->run_name + "_" + sid_._to_string() + "_" +
                            sparams_->name + "." + GetAnalysisName() +
                            ".analysis";
    Logger::Debug("Initializing analysis output file %s", file_name.c_str());
    output_.open(file_name, std::ios::out);
  }
  virtual void InitAnalysis() {}
  virtual void RunAnalysis() {}
  virtual void EndAnalysis() {}

public:
  virtual void Init(std::vector<T> &members, species_parameters<S> &sparams) {
    AnalysisBase::Init();
    sid_ = species_id::_from_integral(S);
    members_ = &members;
    sparams_ = &sparams;
    n_members_ = members_->size();
    InitOutput();
    InitAnalysis();
    Logger::Info("Running %s analysis", GetAnalysisName().c_str());
  }
  virtual void Run() {
    time_ = 0.5 * params_->i_step * params_->delta;
    RunAnalysis();
    iteration_++;
  }
  virtual void End() {
    EndAnalysis();
    AnalysisBase::End();
  }
};

#endif // _SIMCORE_ANALYSIS_H_
