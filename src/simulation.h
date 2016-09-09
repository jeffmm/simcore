#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "space.h"
#include "auxiliary.h"
#include "write_outputs.h"
#include "graphics.h"
#include "species.h"
#include "objects.h"
#include "uberengine.h"
#include "helpers.h"
#include "output_manager.h"

class Simulation {

  private:
    int i_step_,
        n_steps_;
    double time_,
           cpu_init_time_;
    std::string run_name_;
    std::string posit_file_;
    std::ifstream ip;

    OutputManager output_mgr_;
    system_parameters params_;
    rng_properties rng_;
    
    Graphics graphics_;
    SpaceProperties space_;
    UberEngine uengine_;
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
    void IntegrateMP();
    void InteractMP();
    void KineticMonteCarloMP();
    void ZeroForces();
    void DumpAll(int i_step);
    std::vector<graph_struct*> graph_array;

  public:
    Simulation();
    ~Simulation();
    void Run(system_parameters params, std::string name);
    void CreateMovie(system_parameters params, std::string name, std::string posit_file);
};

#endif // _SIMCORE_SIMULATION_H_  
