#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "space.h"
#include "graphics.h"
#include "species.h"
#include "output_manager.h"
#include "interaction_engine.h"
#include "auxiliary.h"
#include "helpers.h"
#include "objects.h"

class Simulation {

  private:
    int i_step_ = 0,
        n_steps_;
    double time_,
           cpu_init_time_;
    std::string run_name_;
    std::vector<std::string> posit_files_;

    OutputManager output_mgr_;
    system_parameters params_;
    rng_properties rng_;

    InteractionEngine iengine_;

    #ifndef NOGRAPH
    Graphics graphics_;
    #endif
    Space space_;
    std::vector<SpeciesBase*> species_;
    rfh::factory species_factory_;
    void InitSimulation();
    void InitSpecies();
    void InitPositInput();
    void ClearSpecies();
    void InitOutputs();
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

  public:
    Simulation() {}
    void Run(system_parameters params, std::string name);
    void CreateMovie(system_parameters params, std::string name, std::vector<std::string> posit_files);
};

#endif // _SIMCORE_SIMULATION_H_  
