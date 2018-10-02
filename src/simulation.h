#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "space.h"
#include "graphics.h"
#include "species.h"
#include "output_manager.h"
#include "interaction_engine.h"
#include "auxiliary.h"
#include "helpers.h"
#include "filament.h"
#include "centrosome.h"
#include "bead_spring.h"
#include "spherocylinder.h"
#include "spindle.h"
#include "passive_filament.h"
//#include "br_bead.h"
//#include "objects.h"
#include "parse_flags.h"

class Simulation {

  private:
    int i_step_ = 0,
        n_steps_,
        frame_num_ = 0;
    double time_,
           cpu_init_time_;
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
    std::vector<SpeciesBase*> species_;
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
    void Draw();
    void WriteOutputs();
    void GetGraphicsStructure();
    void Integrate();
    void Interact();
    void ReadSpeciesPositions();
    void ZeroForces();
    void Statistics();
    void ScaleSpeciesPositions();
    std::vector<graph_struct*> graph_array;
    void PrintComplete();
    void InsertSpecies(bool force_overlap = false, bool processing = false);
    void RunProcessing(run_options run_opts);
    void InitGraphics();
    void InitProcessing(run_options run_opts);

  public:
    Simulation() {}
    void Run(system_parameters params);
    void ProcessOutputs(system_parameters params, run_options run_opts);
};

#endif // _SIMCORE_SIMULATION_H_  
