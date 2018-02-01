
#include "simulation.h"

#define REGISTER_SPECIES(n,m) species_factory_.register_class<n>(#m);

void Simulation::Run(system_parameters params) {
  params_ = params;
  run_name_ = params.run_name;
  InitSimulation();
  RunSimulation();
  ClearSimulation();
}

void Simulation::RunSimulation() {
  std::cout << "  Running simulation" << std::endl;
  for (i_step_ = 0; i_step_<params_.n_steps; ++i_step_) {
    time_ = (i_step_+1) * params_.delta; 
    PrintComplete();
    ZeroForces();
    Interact();
    Integrate();
    Statistics();
    Draw();
    if (early_exit) {
      early_exit = false;
      std::cout << "  Early exit triggered. Ending simulation.\n";
      return;
    }
    WriteOutputs();
  }
}

void Simulation::PrintComplete() {
  int it,steps;
  if (params_.n_steps > 10000) {// && i_step_ >= 100) {
    it = i_step_/100;
    steps = params_.n_steps/10000;
  }
  else {
    it = i_step_*100;
    steps = params_.n_steps;
  }
  if (it % (steps/10) == 0) {
    printf("    %2.1f%% complete\r", (double) it/ steps);
    fflush(stdout);
  }
  if (debug_trace) {
    printf("********\nStep %d\n********\n", i_step_);
  }
}


void Simulation::Integrate() {
  for (auto it=species_.begin(); it!=species_.end(); ++it) {
    (*it)->UpdatePositions();
  }
}

void Simulation::Interact() {
  iengine_.Interact();
}

void Simulation::ZeroForces() {
  for (auto it=species_.begin(); it != species_.end(); ++it) {
    (*it)->ZeroForces();
  }
}

void Simulation::Statistics() {
  if (i_step_ % params_.n_thermo == 0 && i_step_ > 0) {
    iengine_.CalculatePressure();
    if (params_.constant_pressure) {
      space_.ConstantPressure();
    }
    else if (params_.constant_volume) {
      space_.ConstantVolume();
    }
  }
  if (space_.GetUpdate()) {
    space_.UpdateSpace();
    ScaleSpeciesPositions();
  }
}

void Simulation::ScaleSpeciesPositions() {
  for (auto spec : species_) {
    spec->ScalePositions();
  }
}

void Simulation::InitSimulation() {
  std::cout << "  Initializing simulation" << std::endl;
  rng_.Init(params_.seed);
  space_.Init(&params_);
  InitObjects();
  InitSpecies();
  iengine_.Init(&params_, &species_, space_.GetStruct());
  InsertSpecies(params_.load_checkpoint, params_.load_checkpoint);
  InitOutputs();
  if (params_.graph_flag) {
    InitGraphics();
  }
}

void Simulation::InitObjects() {
  Object::SetParams(&params_);
  Object::SetNDim(params_.n_dim);
  Object::SetDelta(params_.delta);
  Object::SetSeed(gsl_rng_get(rng_.r));
  Object::SetSpace(space_.GetStruct());
}

void Simulation::InitGraphics() {
  GetGraphicsStructure();
  double background_color = (params_.graph_background == 0 ? 0.1 : 1);
  #ifndef NOGRAPH
  graphics_.Init(&graph_array, space_.GetStruct(), background_color, params_.draw_boundary);
  graphics_.DrawLoop();
  #endif
  params_.movie_directory.append("/");
  params_.movie_directory.append(params_.run_name);
  // Grab first frame
  #ifndef NOGRAPH
  if (params_.movie_flag) {
    // Record bmp image of frame 
    grabber(graphics_.windx_, graphics_.windy_,
            params_.movie_directory, (int) i_step_/params_.n_graph);
  }
  #endif
}

void Simulation::InitSpecies() {

  // We have to search for the various types of species that we have
  //REGISTER_SPECIES(MDBeadSpecies,md_bead);
  //REGISTER_SPECIES(HardRodSpecies,hard_rod);
  //REGISTER_SPECIES(BrBeadSpecies,br_bead);
  REGISTER_SPECIES(CentrosomeSpecies,centrosome);
  REGISTER_SPECIES(FilamentSpecies,filament);
  REGISTER_SPECIES(BeadSpringSpecies,bead_spring);
  REGISTER_SPECIES(SpherocylinderSpecies,spherocylinder);

  /* Search the species_factory_ for any registered species,
   and find them in the yaml file */
  species_.reserve(species_factory_.m_classes.size());
  for (auto registered = species_factory_.m_classes.begin();
      registered != species_factory_.m_classes.end(); ++registered) {
    SpeciesBase *spec = (SpeciesBase*)species_factory_.construct(registered->first);
    spec->Init(&params_, space_.GetStruct(), gsl_rng_get(rng_.r));
    if (spec->GetNInsert() > 0) {
      species_.push_back(spec);
      species_.back()->Reserve();
    }
  }
}

void Simulation::InsertSpecies(bool force_overlap, bool processing) {
  printf("\r  Inserting species: 0%% complete");
  fflush(stdout);
  // Assuming Random insertion for now
  //force_overlap = true;
  for (auto spec = species_.begin(); spec!=species_.end(); ++spec) {
  // Check for random insertion
    if (processing || (*spec)->GetInsertionType().find("random") == std::string::npos) {
      // Insertion not random, force overlap
      force_overlap = true;
    }
    int num = (*spec)->GetNInsert();
    bool not_done = true;
    int inserted = 0;
    int num_attempts = 0;
    while(num != inserted) {
      inserted = 0;
      int num_failures = 0;
      while(num != inserted) {
        (*spec)->AddMember();
        inserted++;
        // First check that we are respecting boundary conditions
        if (params_.boundary != 0 && !processing && iengine_.CheckBoundaryConditions()) {
          (*spec)->PopMember();
          inserted--;
          // We are not counting boundary condition failures in insertion
          // failures, since insertion failures are for packing issues
        }
        // Check if we have an overlap of objects
        else if (!force_overlap && !(*spec)->CanOverlap() && !processing && iengine_.CheckOverlap()) {
          //printf("Has an overlap\n");
          (*spec)->PopMember();
          inserted--;
          num_failures++;
        }
        // Otherwise update display of percentage of species inserted
        else {
          printf("\r  Inserting species: %d%% complete", (int)(100 * (float)inserted / (float)num));
          fflush(stdout);
        }
        if (num_failures>params_.species_insertion_failure_threshold) {
          break;
        }
      }
      putchar('\n');
      if (num != inserted) {
        printf("  Species insertion failure threshold reached. Reattempting insertion.\n");
        (*spec)->PopAll();
      }
      if (++num_attempts > 20) {
        error_exit("Unable to insert species randomly within the hard-coded attempt threshold of 20!\n");
      }
    }
    if (!processing) {
      printf("\n");
      if ((*spec)->GetInsertionType().find("random") == std::string::npos) {
        (*spec)->ArrangeMembers();
        if (!(*spec)->CanOverlap() && iengine_.CheckOverlap()) {
          error_exit("Species inserted with deterministic insertion type is overlapping!");
        }
      }
    }
  }
}

void Simulation::ClearSpecies() {
  for (auto it=species_.begin(); it!=species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  species_factory_.clear();
}

void Simulation::ClearSimulation() {
  output_mgr_.Close();
  ClearSpecies();
  #ifndef NOGRAPH
  if (params_.graph_flag) {
    graphics_.Clear();
  }
  #endif
  std::cout << "  Simulation complete" << std::endl;
}

void Simulation::Draw() {
  #ifndef NOGRAPH
  if (params_.graph_flag && i_step_%params_.n_graph==0) {
    GetGraphicsStructure();
    graphics_.Draw();
    if (params_.movie_flag) {
      // Record bmp image of frame 
      grabber(graphics_.windx_, graphics_.windy_,
              params_.movie_directory, (int) i_step_/params_.n_graph);
    }
  }
  #endif
}

void Simulation::GetGraphicsStructure() {
  graph_array.clear();
  for (auto it=species_.begin(); it!=species_.end(); ++it) {
    (*it)->Draw(&graph_array);
  }
}

void Simulation::InitOutputs() {
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), &i_step_, run_name_);
  if (params_.time_flag) {
    cpu_init_time_ = cpu_time();
  }
}

void Simulation::InitInputs(bool posits_only, int reduce_factor) {
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), &i_step_, run_name_, true, posits_only, reduce_factor);
}

void Simulation::WriteOutputs() {
  output_mgr_.WriteOutputs();
  if (i_step_ == 0) {
    return; // skip first step
  }
  if (params_.time_flag && i_step_ == params_.n_steps-1) {
    double cpu_final_time = cpu_time();
    double cpu_time = cpu_final_time - cpu_init_time_;
    std::cout << "CPU Time for Initialization: " <<  cpu_init_time_ << "\n";
    std::cout << "CPU Time: " << cpu_time << "\n";
    std::cout << "Sim Time: " << time_ << "\n";
    std::cout << "CPU Time/Sim Time: " << "\n" << cpu_time/time_ << std::endl;
  }
}

void Simulation::ProcessOutputs(system_parameters params, run_options run_opts) {
  params_ = params;
  run_name_ = params.run_name;
  InitProcessing(run_opts);
  RunProcessing(run_opts.analysis_flag);
  ClearSimulation();
}

// Initialize everything we need for processing
void Simulation::InitProcessing(run_options run_opts) {
  rng_.Init(params_.seed);
  space_.Init(&params_);
  InitObjects();
  InitSpecies();
  InsertSpecies(true, true);
  if (run_opts.reduce_flag) {
    InitInputs(run_opts.use_posits,run_opts.reduce_factor);
  }
  else {
    InitInputs(run_opts.use_posits, 1);
  }
  if (run_opts.graphics_flag || run_opts.make_movie) {
    params_.graph_flag = 1;
    if (run_opts.use_posits && params_.n_graph < output_mgr_.GetNPosit()) {
      params_.n_graph = output_mgr_.GetNPosit();
    }
    else if (!run_opts.use_posits && params_.n_graph < output_mgr_.GetNSpec()) {
      params_.n_graph = output_mgr_.GetNSpec();
    }
    if (run_opts.make_movie) {
      params_.movie_flag = 1;
    }
    InitGraphics();
  }
  else {
    params_.graph_flag = 0;
  }
  if (run_opts.analysis_flag) {
    for (auto it=species_.begin(); it!=species_.end(); ++it) {
      (*it)->InitAnalysis();
    }
  }
}

void Simulation::RunProcessing(int run_analyses) {
  std::cout << "Processing outputs for: " << run_name_ << std::endl;
  for (i_step_ = 1; i_step_<params_.n_steps; ++i_step_) {
    time_ = (i_step_+1) * params_.delta; 
    PrintComplete();
    output_mgr_.ReadInputs(); 
    if (early_exit) {
      early_exit = false;
      std::cout << "  Early exit triggered. Ending simulation.\n";
      return;
    }
    Draw();
    if (run_analyses) {
      for (auto it=species_.begin(); it!=species_.end(); ++it) {
        if ( ((*it)->GetPositFlag() && i_step_%(*it)->GetNPosit()==0) 
            || ((*it)->GetSpecFlag() && i_step_%(*it)->GetNSpec()==0) ) {
          (*it)->RunAnalysis();
        }
      }
    }
  }
  if (run_analyses) {
    for (auto it=species_.begin(); it!=species_.end(); ++it) {
      (*it)->FinalizeAnalysis();
    }
  }
  // ClearSimulation will run CloseFiles which runs
  // FinalizeAnalysis for each species
}

