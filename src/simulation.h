#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "anchor_list_generic.h"
#include "space.h"
#include "auxiliary.h"
#include "write_outputs.h"
#include "graphics.h"
#include "species.h"
#include "objects.h"
#include "uberengine.h"
#include "uberengine_v2.h"
#include "helpers.h"
#include "output_manager.h"

class Simulation {

  private:
    int i_step_,
        n_steps_;
    int utype_;
    double time_,
           cpu_init_time_;
    std::string run_name_;
    std::vector<std::string> posit_files_;
    //std::fstream ip_;

    OutputManager output_mgr_;
    system_parameters params_;
    rng_properties rng_;
   
    #ifndef NOGRAPH
    Graphics graphics_;
    #endif
    SpaceProperties space_;
    UberEngine uengine_;
    UberEngineV2 uenginev2_;
    std::vector<SpeciesBase*> species_;
    rfh::factory species_factory_;
    al_set anchors_;
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
    void ReadSpeciesPositions();
    void KineticMonteCarloMP();
    void ZeroForces();
    void DumpAll(int i_step);
    std::vector<graph_struct*> graph_array;
    Simulation(const Simulation&){};
    void operator= (Simulation&);

    void ConfigureSpindle();

  public:
    Simulation();
    ~Simulation();
    void Run(system_parameters params, std::string name);
    void CreateMovie(system_parameters params, std::string name, std::vector<std::string> posit_files);
};

#endif // _SIMCORE_SIMULATION_H_  
