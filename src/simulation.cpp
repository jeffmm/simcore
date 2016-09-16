
#include "simulation.h"

#define REGISTER_SPECIES(n,m) species_factory_.register_class<n>(#m);

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
  std::cout << "    steps: " << params_.n_steps << std::endl;
  for (i_step_=0; i_step_<params_.n_steps; ++i_step_) {
    time_ = (i_step_+1) * params_.delta; 
    //if (i_step_ % (params_.n_steps / 100) == 0) {
    if ((100*i_step_) % (params_.n_steps) == 0) {
      printf("%d%% Complete\n", (int)(100 * (float)i_step_ / (float)params_.n_steps));
      fflush(stdout);
    }
    if (debug_trace)
      printf("********\nStep %d\n********\n", i_step_);
    ZeroForces();
    KineticMonteCarloMP();
    InteractMP();
    IntegrateMP();
    // Only will run if DEBUG is enabled
    #ifdef DEBUG
    if (debug_trace)
      DumpAll(i_step_);
    #endif
    Draw();
    WriteOutputs();
  }
}

void Simulation::RunMovie(){
  std::cout << "Running movie: " << run_name_ << "\n";
  std::cout << "    steps: " << params_.n_steps << std::endl;
  for (i_step_=0; i_step_<params_.n_steps; ++i_step_) {
    time_ = (i_step_+1) * params_.delta; 
    //if (i_step_ % (params_.n_steps / 100) == 0) {
    if ((100*i_step_) % (params_.n_steps) == 0) {
      printf("%d%% Complete\n", (int)(100 * (float)i_step_ / (float)params_.n_steps));
      fflush(stdout);
    }
    if (debug_trace)
      printf("********\nStep %d\n********\n", i_step_);
    if (i_step_%params_.n_posit == 0)
      ReadSpeciesPositions(); 
    
    //ZeroForces();
    //KineticMonteCarloMP();
    //InteractMP();
    //IntegrateMP();
    // Only will run if DEBUG is enabled
    //#ifdef DEBUG
    //if (debug_trace)
      //DumpAll(i_step_);
    //#endif
    Draw();
    //WriteOutputs();
  }
}

void Simulation::DumpAll(int i_step) {
    // Very yucky dump of all the particles and their positions and forces
    uengine_.DumpAll();
}

void Simulation::Integrate() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->UpdatePositions();
}

void Simulation::IntegrateMP() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->UpdatePositionsMP();
}

void Simulation::InteractMP() {
  uengine_.InteractMP();
}

void Simulation::KineticMonteCarloMP() {
  uengine_.StepKMC();
}

void Simulation::ZeroForces() {
  for (auto it=species_.begin(); it != species_.end(); ++it) {
    (*it)->ZeroForces();
  }
}

void Simulation::InitSimulation() {

  space_.Init(&params_, gsl_rng_get(rng_.r));
  output_mgr_.Init(&params_, &i_step_);
  InitSpecies();
  uengine_.Init(&params_, space_.GetStruct(), &species_, gsl_rng_get(rng_.r));
  if (params_.graph_flag) {
    GetGraphicsStructure();
    double background_color = (params_.graph_background == 0 ? 0.1 : 1);
    graphics_.Init(&graph_array, space_.GetStruct(), background_color);
    graphics_.DrawLoop();
  }
  InitOutputs();
}


void Simulation::InitSpecies() {
//#include "init_species.h"
  // Check out the configuration file
  std::cout << "********\n";
  std::cout << "Species Load ->\n";
  std::cout << "  file: " << params_.config_file << std::endl;

  YAML::Node node = YAML::LoadFile(params_.config_file);

  // We have to search for the various types of species that we have
  // Maybe hijack init_species.h for this
  REGISTER_SPECIES(MDBeadSpecies,md_bead);
  REGISTER_SPECIES(BrRodSpecies,br_rod);
  REGISTER_SPECIES(XlinkSpecies,xlink);
  REGISTER_SPECIES(FilamentSpecies,filament);
  REGISTER_SPECIES(MDBeadOptSpecies,md_bead_opt);
  REGISTER_SPECIES(BrBeadSpecies,br_bead);

  // Search the species_factory_ for any registered species, and find them in the
  // yaml file
  for (auto possibles = species_factory_.m_classes.begin(); possibles != species_factory_.m_classes.end(); ++possibles) {
    if (node[possibles->first]) {
      SpeciesBase *spec = (SpeciesBase*)species_factory_.construct(possibles->first);
      spec->InitConfig(&params_, space_.GetStruct(), gsl_rng_get(rng_.r), &output_mgr_);
      spec->Configurator();
      species_.push_back(spec);
      output_mgr_.AddSpecie(spec);
    }
  }
}

//FIXME Only works for one species at a time -AL
  //TODO Have output_manager run movies
void Simulation::InitPositInput(){
  for (auto pos_it : posit_files_){
    int nchar;
    std::fstream ip;

    ip.open(pos_it, std::ios::binary | std::ios::in );
    if (!ip.is_open()){
      std::cout<<"Input "<< pos_it <<" file did not open\n";
      exit(1);
    }
    else{
      //Get Sim data from posit file
      ip.read(reinterpret_cast<char*>(&nchar), sizeof(int));
      std::string sid_str(nchar, ' ');
      ip.read(&sid_str[0], nchar);
      SID sid_posit = StringToSID(sid_str);
      ip.read(reinterpret_cast<char*>(&(params_.n_steps)), sizeof(int));
      ip.read(reinterpret_cast<char*>(&(params_.n_posit)), sizeof(int));

      //Get the location where species data starts and close input for now
      std::ios::streampos beg = ip.tellg();
      ip.close();

      for (auto spec_it : species_ )
        if (sid_posit == spec_it->GetSID()){
          std::cout<<"SID of posit file is "<< sid_str << std::endl;
          spec_it->InitInputFile(pos_it, beg);
          break;
        }
    }
  }
}

void Simulation::ClearSpecies() {
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    delete (*it);
}

void Simulation::ClearSimulation() {
  space_.Clear();
  if (params_.posit_flag == 1)
    output_mgr_.Close();
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
    uengine_.Draw(&graph_array);
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

  {
    uengine_.PrepOutputs();
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

  output_mgr_.WriteOutputs();

  if (i_step_ == params_.n_steps-1) {
    for (auto it=species_.begin(); it!=species_.end(); ++it)
      (*it)->WriteOutputs(run_name_);
  }

  // XXX CJE FIXME write outputs more clearly
  if (i_step_%1000==0) {
    uengine_.WriteOutputs(i_step_);
  }

}

void Simulation::CreateMovie(system_parameters params, std::string name, std::vector<std::string> posit_files){
  params_ = params;

  //Graph and don't make new posit files
  params_.graph_flag = 1;
  params_.posit_flag = 0;

  run_name_ = name;
  posit_files_ = posit_files;
  rng_.init(params_.seed);
  InitSimulation();
  InitPositInput();
  RunMovie();
  ClearSimulation();
}

void Simulation::ReadSpeciesPositions(){
  for (auto it=species_.begin(); it!=species_.end(); ++it)
    (*it)->ReadPosits();
}
