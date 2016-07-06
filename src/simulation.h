#ifndef _SIMCORE_SIMULATION_H_
#define _SIMCORE_SIMULATION_H_

#include "space.h"
#include "auxiliary.h"
#include "forces.h"
#include "write_outputs.h"
#include "graphics.h"
#include "species.h"
#include "objects.h"

class Simulation {

  private:
    int i_step_;
    double time_,
           cpu_init_time_;
    std::string run_name_;
    system_parameters params_;
    rng_properties rng_;
    
    Graphics graphics_;
    SpaceProperties space_;
    Forces forces_;
    std::vector<SpeciesBase*> species_;
    void InitSimulation();
    void InitSpecies();
    void ClearSpecies();
    void InitOutputs();
    void RunSimulation();
    void ClearSimulation();
    void Draw();
    void WriteOutputs();
    void GetGraphicsStructure();
    void Integrate();
    void IntegrateMP();
    void InteractMP();
    void KineticMonteCarloMP();
    void DumpAll(int i_step);
    std::vector<graph_struct*> graph_array;

  public:
    Simulation();
    ~Simulation();
    void Run(system_parameters params, std::string name);
};

#endif // _SIMCORE_SIMULATION_H_  
