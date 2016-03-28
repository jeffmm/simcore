#ifndef _CYTOSCORE_SIMULATION_H_
#define _CYTOSCORE_SIMULATION_H_

#include "space.h"
#include "auxiliary.h"
//#include "filament.h"
//#include "sphere.h"
//#include "external_forces.h"
//#include "parameters.h"
//#include "write_outputs.h"
#include "graphics.h"
#include "particle_md_system.h"

class Simulation {

  private:
    int i_step,
        i_run;
    double time,
           cpu_init_time;
    system_parameters *params;
    rng_properties rng;
    
    Graphics graphics;
    graph_struct g_struct;
    SpaceProperties space;
    //IntegratorManager integrator;
    //ExternalForces forces;
    std::vector<ObjectSystemBase*> systems_;
    void InitSimulation();
    void InitSpace();
    void InitSystems();
    void InitGraphics();
    void InitOutputs();
    void RunSimulation();
    void ClearSimulation();
    void ClearSpace();
    void ClearGraphics();
    void Draw();
    void WriteOutputs();
    void GetGraphicsStructure();

  public:
    Simulation(); // Debug mode
    Simulation(std::string param_file, std::string run_name, int n_runs); // Output files with prefix run_name
};

#endif // _CYTOSCORE_SIMULATION_H_  
