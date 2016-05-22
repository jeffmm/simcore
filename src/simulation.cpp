
#include "simulation.h"

unsigned int Simple::next_oid_ = 0;

Simulation::Simulation() {}
Simulation::~Simulation() {}

void Simulation::Run(system_parameters params, std::string name) {
  params_ = params;
  run_name_ = name;
  rng_.init(params_.seed);
  InitSimulation();
  RunSimulation();
  ClearSimulation();
}

void Simulation::RunSimulation() {
  std::cout << "Running simulation: " << run_name_ << "\n";

  for (i_step_=0; i_step_<params_.n_steps; ++i_step_) {
    time_ = (i_step_+1) * params_.delta; 
    Integrate();
    Draw();
    WriteOutputs();
  }
}

void Simulation::Integrate() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->UpdatePositions();
}

void Simulation::InitSimulation() {

  space_.Init(&params_, gsl_rng_get(rng_.r));
  InitSpecies();
  if (params_.graph_flag) {
    GetGraphicsStructure();
    double background_color = (params_.graph_background == 0 ? 0.1 : 1);
    graphics_.Init(&graph_array, space_.GetStruct(), background_color);
    graphics_.DrawLoop();
  }
}

void Simulation::InitSpecies() {
#include "init_species.cpp"
}
void Simulation::ClearSpecies() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    delete (*it);
}

void Simulation::ClearSimulation() {
  space_.Clear();
  ClearSpecies();
  if (params_.graph_flag)
    graphics_.Clear();
}

void Simulation::Draw() {
  if (params_.graph_flag && i_step_%params_.n_graph==0) {
    GetGraphicsStructure();
    graphics_.Draw();
    // Record bmp image of frame 
    if (params_.grab_flag) {
      grabber(graphics_.windx_, graphics_.windy_,
              params_.grab_file, (int) i_step_/params_.n_graph);
    }
  }
}

void Simulation::GetGraphicsStructure() {

  graph_array.clear();
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->Draw(&graph_array);
}

void Simulation::InitOutputs() {
  if (params_.time_flag) {
    cpu_init_time_ = cpu();
  }
}

void Simulation::WriteOutputs() {
  if (i_step_==0) {
    InitOutputs();
    return;
  }
  if (params_.time_flag && i_step_ == params_.n_steps-1) {
    double cpu_time = cpu() - cpu_init_time_;
    std::cout << "CPU Time for Initialization: " <<  cpu_init_time_ << std::endl;
    std::cout << "CPU Time: " << cpu_time << std::endl;
    std::cout << "Sim Time: " << time_ << std::endl;
    std::cout << "CPU Time/Sim Time: " << std::endl << cpu_time/time_ << std::endl;
  }
}

