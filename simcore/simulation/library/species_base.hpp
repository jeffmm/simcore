#ifndef _SIMCORE_SPECIES_BASE_H_
#define _SIMCORE_SPECIES_BASE_H_

#include "auxiliary.hpp"
#include "object.hpp"
#include "yaml-cpp/yaml.h"

class SpeciesBase {
private:
  species_id sid_;

protected:
  int n_members_ = 0, spec_file_iterator_ = -1;
  std::string spec_name_;
  system_parameters *params_;
  space_struct *space_;
  RNG rng_;
  std::fstream oposit_file_;
  std::fstream iposit_file_;
  std::fstream ospec_file_;
  std::fstream ispec_file_;
  std::string checkpoint_file_;
  std::string spec_label_ = "generic";
  std::vector<std::string> spec_file_names_;

public:
  SpeciesBase() {}
  void SetSID(species_id sid) {
    sid_ = sid;
    spec_name_ = sid_._to_string();
  }
  virtual void UpdatePositions() {}
  virtual void Draw(std::vector<graph_struct *> *graph_array) {}
  virtual void Init(system_parameters *params, spec_params *sparams, space_struct *space);
  virtual void ZeroForces() {}
  virtual void GetInteractors(std::vector<Object *> *ix) {}
  virtual void GetLastInteractors(std::vector<Object *> *ix) {}
  virtual double GetPotentialEnergy() { return 0; }
  virtual void ScalePositions() {}
  virtual void AddMember() {}
  virtual void SetLastMemberPosition(double const *const pos) {}
  virtual void PopMember() {}
  virtual void PopAll() {}
  virtual double GetSpecLength() { return 0; }
  virtual double GetSpecDiameter() { return 0; }
  virtual void ArrangeMembers() {}
  virtual int CanOverlap() { return -1; }
  species_id const GetSID() { return sid_; }
  virtual void Report() {}
  int const GetNMembers() { return n_members_; }
  int const GetNInsert() { return -1; }
  int const GetNPosit() { return -1; }
  int const GetNSpec() { return -1; }
  int const GetNCheckpoint() { return -1; }
  bool const GetPositFlag() { return false; }
  bool const GetSpecFlag() { return false; }
  bool const GetCheckpointFlag() { return false; }
  std::string GetInsertionType() { return ""; }
  virtual int GetCount() { return 0; }
  virtual void WriteOutputs(std::string run_name) {}
  virtual void WritePosits() {}
  virtual void WriteSpecs() {}
  virtual void WriteCheckpoints() {}
  virtual void ReadSpecs() {}
  virtual void ReadCheckpoints() {}
  virtual void ReadPosits() {}
  virtual void ReadPositsFromSpecs() {}
  virtual void InitAnalysis() {}
  virtual void RunAnalysis() {}
  virtual void FinalizeAnalysis() {}
  virtual void InitOutputFiles(std::string run_name);
  virtual void InitPositFile(std::string run_name);
  virtual void InitSpecFile(std::string run_name);
  virtual void InitPositFileInput(std::string run_name);
  virtual void InitSpecFileInput(std::string run_name);
  virtual bool InitSpecFileInputFromFile(std::string run_name);
  virtual bool HandleEOF();
  virtual void InitInputFiles(std::string run_name, bool posits_only,
                              bool with_reloads);
  virtual void InitCheckpoints(std::string run_name);
  virtual void LoadFromCheckpoints(std::string run_name,
                                   std::string checkpoint_run_name);
  virtual int OutputIsOpen() { return oposit_file_.is_open(); }
  virtual int InputIsOpen() { return iposit_file_.is_open(); }
  virtual void CloseFiles();
  virtual void CleanUp() {}
  virtual void Reserve() {}
  virtual double const GetVolume() { return 0; }
  virtual double const GetDrMax() { return 0; }
  virtual void ZeroDrTot() {}
  virtual void CustomInsert() {}
  virtual const bool CheckInteractorUpdate() { return false; }
};

#endif
