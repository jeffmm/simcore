#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "auxiliary.hpp"
#include "graphics.hpp"
#include "helpers.hpp"
#include "interaction_engine.hpp"
#include "output_manager.hpp"
#include "space.hpp"
#include "filament_species.hpp"
//#include "centrosome.hpp"
//#include "bead_spring.hpp"
#include "spherocylinder.hpp"
//#include "spindle.hpp"
#include "br_bead.hpp"
//#include "objects.hpp"
#include "parse_flags.hpp"

class Simulation {
private:
  int i_step_ = 0;
  int n_steps_;
  int frame_num_ = 0;
  double time_;
  double cpu_init_time_;
  std::string run_name_;
  std::vector<std::string> posit_files_;

  OutputManager output_mgr_;
  system_parameters params_;
  RNG rng_;

  InteractionEngine iengine_;

#ifndef NOGRAPH
  Graphics graphics_;
#endif
  Space space_;
  std::vector<SpeciesBase *> species_;
  rfh::factory species_factory_;
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
  UNIT_TESTER;

public:
  Simulation() {}
  void Run(system_parameters params);
  void ProcessOutputs(system_parameters params, run_options run_opts);
};

#endif // _SIMCORE_SIMULATION_H_
