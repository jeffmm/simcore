#ifndef _CYTOSCORE_SIMULATION_H_
#define _CYTOSCORE_SIMULATION_H_

#include "space.h"
#include "auxiliary.h"
//#include "filament.h"
//#include "sphere.h"
//#include "external_forces.h"
//#include "parameters.h"
#include "write_outputs.h"
#include "graphics.h"
#include "objects.h"
//#include "particle_md_system.h"

class Simulation {

  private:
    int i_step,
        i_run;
    double time,
           cpu_init_time;
    std::string run_name_;
    system_parameters params_;
    rng_properties rng;
    
    Graphics graphics;
    graph_struct g_struct;
    SpaceProperties space;
    //IntegratorManager integrator;
    //ExternalForces forces;
    //std::vector<ObjectSystemBase*> systems_;
    std::vector<BrownianDimer> dimers_;
    void InitSimulation();
    //void InitSpace();
    void InitSystems();
    //void InitGraphics();
    //void InitOutputs();
    void RunSimulation();
    void ClearSimulation();
    //void ClearSpace();
    //void ClearGraphics();
    void Draw();
    //void WriteOutputs();
    void GetGraphicsStructure();
    void Integrate();
    std::vector<graph_struct*> graph_array;

  public:
    Simulation();
    ~Simulation();
    void Run(system_parameters params, std::string name);
};

#endif // _CYTOSCORE_SIMULATION_H_  
