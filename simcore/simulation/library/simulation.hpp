#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "graphics.hpp"
#include "interaction_manager.hpp"
#include "output_manager.hpp"
#include "space.hpp"

class Simulation {
 private:
  UNIT_TEST
  int i_step_ = 0;
  int log_interval_ = 10;
  int n_steps_;
  int frame_num_ = 0;
  double time_;
  double cpu_init_time_;
  std::string run_name_;
  std::vector<std::string> posit_files_;

  OutputManager output_mgr_;
  system_parameters params_;
  RNG *rng_;
  ParamsParser parser_;

  InteractionManager ix_mgr_;

#ifndef NOGRAPH
  Graphics graphics_;
#endif
  Space space_;
  std::vector<SpeciesBase *> species_;
  void InitSimulation();
  void InitObjects();
  void InitSpecies();
  void InitPositInput();
  void ClearSpecies();
  void InitOutputs();
  void InitInputs(run_options run_opts);
  void RunSimulation();
  void RunMovie();
  void ClearSimulation();
  void Draw(bool single_frame = false);
  void WriteOutputs();
  void GetGraphicsStructure();
  void Integrate();
  void Interact();
  void ReadSpeciesPositions();
  void ZeroForces();
  void Statistics();
  void ScaleSpeciesPositions();
  std::vector<graph_struct *> graph_array_;
  void PrintComplete();
  void InsertSpecies(bool force_overlap = false, bool processing = false);
  void RunProcessing(run_options run_opts);
  void InitGraphics();
  void InitProcessing(run_options run_opts);

 public:
  Simulation() {}
  void Run(YAML::Node sim_params);
  void ProcessOutputs(YAML::Node sim_params, run_options run_opts);
};

#endif  // _SIMCORE_SIMULATION_H_
