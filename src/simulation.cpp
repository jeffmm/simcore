
#include "simulation.h"

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
    //DPRINTF("********\nStep %d\n********\n", i_step_);
    //Interact();
    //Integrate();
    InteractMP();
    IntegrateMP();
    // Only will run if DEBUG is enabled
    #ifdef DEBUG
    DumpAll(i_step_);
    #endif
    Draw();
    WriteOutputs();
  }
}

void Simulation::DumpAll(int i_step) {
    // Very yucky dump of all the particles and their positions and forces
    forces_.DumpAll();
}

void Simulation::Integrate() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->UpdatePositions();
}

void Simulation::IntegrateMP() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->UpdatePositionsMP();
}

void Simulation::Interact() {
  if (i_step_%params_.n_update_cells==0) {
    forces_.UpdateCellList(species_);
  }
  forces_.Interact();
}

void Simulation::InteractMP() {
  if (i_step_ % params_.n_update_cells == 0) {
    forces_.UpdateScheme(species_);
  }
  forces_.InteractMP();
}

void Simulation::InitSimulation() {

  space_.Init(&params_, gsl_rng_get(rng_.r));
  InitSpecies();
  forces_.Init(space_.GetStruct(), species_, params_.ftype, params_.cell_length, params_.draw_interactions);
  if (params_.graph_flag) {
    GetGraphicsStructure();
    double background_color = (params_.graph_background == 0 ? 0.1 : 1);
    graphics_.Init(&graph_array, space_.GetStruct(), background_color);
    graphics_.DrawLoop();
  }
  InitOutputs();
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
    if (params_.grab_flag) {
      // Record bmp image of frame 
      grabber(graphics_.windx_, graphics_.windy_,
              params_.grab_file, (int) i_step_/params_.n_graph);
    }
  }
}

void Simulation::GetGraphicsStructure() {

  graph_array.clear();
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->Draw(&graph_array);
  if (params_.draw_interactions)
    forces_.Draw(&graph_array);
}

void Simulation::InitOutputs() {
  //double tot_en=0;
  //for (auto it=species_.begin(); it!=species_.end(); ++it)
    //tot_en += (*it)->GetTotalEnergy();
  //std::cout << "Initial system energy: " << tot_en << std::endl;
  if (params_.time_flag) {
    cpu_init_time_ = cpu();
  }
  if (params_.energy_analysis_flag) {
    std::ostringstream file_name;
    file_name << run_name_ << "-energy.log";
    std::ofstream en_file(file_name.str().c_str(), std::ios_base::out);
    en_file << "#kinetic  #potential  #total\n";
    en_file.close();
  }

}

void Simulation::WriteOutputs() {
  if (i_step_ == 0) {
    return; // skip first step
  }
  if (params_.time_flag && i_step_ == params_.n_steps-1) {
    double cpu_time = cpu() - cpu_init_time_;
    std::cout << "CPU Time for Initialization: " <<  cpu_init_time_ << "\n";
    std::cout << "CPU Time: " << cpu_time << "\n";
    std::cout << "Sim Time: " << time_ << "\n";
    std::cout << "CPU Time/Sim Time: " << "\n" << cpu_time/time_ << std::endl;
    //double tot_en = 0;
    //for (auto it=species_.begin(); it!=species_.end(); ++it)
      //tot_en += (*it)->GetTotalEnergy();
    //std::cout << "Final system energy: " << tot_en << std::endl;
  }

  if (i_step_%1000==0 && params_.energy_analysis_flag) {
    std::ostringstream file_name;
    file_name << run_name_ << "-energy.log";
    std::ofstream en_file(file_name.str().c_str(), std::ios_base::out | std::ios_base::app);
    en_file.precision(16);
    en_file.setf(std::ios::fixed, std::ios::floatfield);
    double tot_en=0;
    double k_en=0;
    double p_en=0;
    for (auto it=species_.begin(); it!=species_.end(); ++it) {
      k_en += (*it)->GetKineticEnergy();
      p_en += (*it)->GetPotentialEnergy();
      tot_en += (*it)->GetTotalEnergy();
    }
    en_file << k_en << " " << p_en << " " << tot_en << "\n";
    en_file.close();
  }


}

