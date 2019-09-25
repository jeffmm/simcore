
#include "simulation.hpp"

#define REGISTER_SPECIES(n, m) species_factory_.register_class<n>(#m);

/* Initialize simulation parameters and run simulation */
void Simulation::Run(system_parameters params) {
  params_ = params;
  run_name_ = params.run_name;
  // Initialize simulation data structures
  InitSimulation();
  // Begin simulation
  RunSimulation();
  // Tear down data structures (e.g. cell lists, etc)
  ClearSimulation();
}

/* Continuously loop over object data structures, solving their equations of
 * motion using Runge-Kutta-like integration schemes for n_steps time steps. */
void Simulation::RunSimulation() {
#ifdef ENABLE_OPENMP
  std::cout << "  Running simulation on " << omp_get_max_threads() << " threads"
            << std::endl;
#else
  std::cout << "  Running simulation" << std::endl;
#endif
  /* Solve equations of motion for all objects */

  /* Note that we are beginning with i_step = 1 here as a convention. This is
   because filaments use a midstep method, so that at the end of step 1 and
   every odd step, they will be in a midstep configuration, whereas at the end
   of every even step, they will be in the fullstep configuration. Since most
   sane step frequency functions (e.g. n_graph, n_thermo, etc) are even numbers,
   we want these to represent filaments in their fullstep configurations. We
   also update other objects' positions on every even step too (e.g. xlinks). */
  for (i_step_ = 1; i_step_ < params_.n_steps + 1; ++i_step_) {
    params_.i_step = i_step_;
    time_ = (i_step_ + 1) * params_.delta;
    // Output progress
    PrintComplete();
    /* Zero the force_ array on all objects for bookkeeping */
    ZeroForces();
    /* Calculate forces between objects in system */
    Interact();
    /* Integrate EOM on all objects to update positions, etc */
    Integrate();
    /* Update system pressure, volume, etc if necessary */
    Statistics();
    /* Draw the particles in graphics window, if using graphics */
    Draw();
    /* Catch soft exceptions and terminate the simulation early */
    if (early_exit) {
      early_exit = false;
      std::cout << "  Early exit triggered. Ending simulation.\n";
      return;
    }
    /* Generate all output files */
    WriteOutputs();
  }
}

/* Print simulation progress as a percentage. If the print_complete parameter
 * evaluates to true, this will print a new line for each progress update,
 * which can be useful for logging. Otherwise, the progress will be refreshed
 * in a single line printed to standard output. */
void Simulation::PrintComplete() {
  int it, steps;
  // Avoid weird issues with longs
  if (params_.n_steps > 100000) {
    it = i_step_ / 100;
    steps = params_.n_steps / 10000;
  } else {
    it = i_step_ * 100;
    steps = params_.n_steps;
  }
  if (params_.print_complete && it % steps == 0) {
    printf("    %d%% complete\n", it / steps);
  } else if (!params_.print_complete && it % steps == 0) {
    printf("    %2.1f%% complete\r", (double)it / steps);
    fflush(stdout);
  }
  if (debug_trace) {
    printf("********\nStep %d\n********\n", i_step_);
  }
}

/* Update the positions of all objects in the system using numerical
 * integration for one time step defined by the delta parameter. */
void Simulation::Integrate() {
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->UpdatePositions();
  }
}

/* Calculate interaction forces between all objects if necessary. */
void Simulation::Interact() { iengine_.Interact(); }

/* Remove forces on all objects. Objects with inertia that use higher order RK
 * schemes with dependence on previous forces need to save their forces in
 * prev_forces_ arrays. */
void Simulation::ZeroForces() {
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->ZeroForces();
  }
}

/* Update system pressure, volume and rescale system size if necessary,
 * handling periodic boundaries in a sane way. */
void Simulation::Statistics() {
  if (i_step_ % params_.n_thermo == 0 && i_step_ > 0) {
    /* Calculate system pressure from stress tensor */
    iengine_.CalculatePressure();
    if (params_.constant_pressure) {
      space_.ConstantPressure();
    }
    /* Constant volume is used to adiabatically adjust the system size to
     achieve high densities from random insertions. */
    else if (params_.constant_volume) {
      space_.ConstantVolume();
    }
  }
  /* Update object positions, etc and handle periodic BCs */
  if (space_.GetUpdate()) {
    space_.UpdateSpace();
    ScaleSpeciesPositions();
  }
}

/* Update object positions in periodic space if the system size changed. */
void Simulation::ScaleSpeciesPositions() {
  for (auto spec : species_) {
    spec->ScalePositions();
  }
}

/* Initialize all the data structures in the simulation */
void Simulation::InitSimulation() {
  std::cout << "  Initializing simulation" << std::endl;
  space_.Init(&params_);
  InitObjects();
  InitSpecies();
  iengine_.Init(&params_, &species_, space_.GetStruct(), &i_step_);
  InsertSpecies(params_.load_checkpoint, params_.load_checkpoint);
  InitOutputs();
  if (params_.graph_flag) {
    InitGraphics();
  }
}

/* Initialize static object parameters that are used everywhere */
void Simulation::InitObjects() {
  Object::SetParams(&params_);
  Object::SetNDim(params_.n_dim);
  Object::SetDelta(params_.delta);
  Object::SetSpace(space_.GetStruct());
}

/* Generate graphics window and draw initial simulation setup */
void Simulation::InitGraphics() {
  GetGraphicsStructure();
  double background_color = (params_.graph_background == 0 ? 0.1 : 1);
// If NOGRAPH is defined, skip drawing and grabbing
#ifndef NOGRAPH
  // Initialize graphics structures
  graphics_.Init(&graph_array_, space_.GetStruct(), background_color,
                 params_.draw_boundary, params_.auto_graph);
  graphics_.DrawLoop();
#endif
  // Initialize directory for grabbed images
  params_.movie_directory.append("/");
  params_.movie_directory.append(params_.run_name);
#ifndef NOGRAPH
  // Grab first frame
  if (params_.movie_flag) {
    // Record bmp image of frame into movie_directory
    grabber(graphics_.windx_, graphics_.windy_, params_.movie_directory,
            (int)i_step_ / params_.n_graph);
  }
#endif
}

/* Initialize object types */
void Simulation::InitSpecies() {
  /* Factories for creating and initializing registered species with their
     assigned species id */
  REGISTER_SPECIES(FilamentSpecies, filament);
  REGISTER_SPECIES(PassiveFilamentSpecies, passive_filament);
  REGISTER_SPECIES(BeadSpringSpecies, bead_spring);
  REGISTER_SPECIES(SpherocylinderSpecies, spherocylinder);
  REGISTER_SPECIES(BrBeadSpecies, br_bead);

  /* Search the species_factory_ for any registered species,
   and find them in the yaml file */
  species_.reserve(species_factory_.m_classes.size());
  for (auto registered = species_factory_.m_classes.begin();
       registered != species_factory_.m_classes.end(); ++registered) {
    SpeciesBase *spec =
        (SpeciesBase *)species_factory_.construct(registered->first);
    spec->Init(&params_, space_.GetStruct(), gsl_rng_get(rng_.r));
    if (spec->GetNInsert() > 0) {
      species_.push_back(spec);
      species_.back()->Reserve();
    }
  }
}

/* Initialize object positions and orientations.*/
void Simulation::InsertSpecies(bool force_overlap, bool processing) {
  if (params_.print_complete) {
    printf("  Inserting species: 0%% complete\n");
  } else {
    printf("\r  Inserting species: 0%% complete");
    fflush(stdout);
  }
  // Assuming Random insertion for now
  // force_overlap = true;
  for (auto spec = species_.begin(); spec != species_.end(); ++spec) {
    // Check for random insertion
    if (processing ||
        (*spec)->GetInsertionType().find("random") == std::string::npos) {
      // Insertion not random, force overlap
      force_overlap = true;
    }
    int num = (*spec)->GetNInsert();
    // printf("n_insert: %d\n",num);
    // exit(0);
    bool not_done = true;
    int inserted = 0;
    int num_attempts = 0;
    while (num != inserted) {
      inserted = 0;
      int num_failures = 0;
      while (num != inserted) {
        (*spec)->AddMember();
        // First check that we are respecting boundary conditions
        std::vector<Object *> last_ixors;
        (*spec)->GetLastInteractors(&last_ixors);
        if (params_.boundary != 0 && !processing &&
            iengine_.CheckBoundaryConditions(last_ixors)) {
          (*spec)->PopMember();
          /* We are not counting boundary condition failures in insertion
           failures, since insertion failures are for packing issues */
          // num_failures++;
        }
        // Check if we have an overlap of objects
        else if (!force_overlap && !(*spec)->CanOverlap() && !processing &&
                 iengine_.CheckOverlap(last_ixors)) {
          (*spec)->PopMember();
          num_failures++;
        }
        // Otherwise update display of percentage of species inserted
        else {
          inserted++;
          if (!force_overlap && !processing) {
            iengine_.AddInteractors(last_ixors);
          }
          int insert_percentage = (int)(100 * (float)inserted / (float)num);
          if (params_.print_complete) {
            if (insert_percentage % 10 == 0) {
              printf("  Inserting species: %d%% complete\n", insert_percentage);
            }
          } else {
            printf("\r  Inserting species: %d%% complete", insert_percentage);
            fflush(stdout);
          }
        }
        if (num_failures > params_.species_insertion_failure_threshold) {
          break;
        }
        /* Update the number of particles we need to insert, in case a species
           needs to have a certain packing fraction */
        num = (*spec)->GetNInsert();
      }
      if (num != inserted) {
        // Attempt a lattice-based insertion strategy (only 2d for now)
        if (params_.n_dim == 3)
          continue;
        warning("Attempting lattice-based insertion strategy");
        double pos[3] = {0, 0, 0};
        pos[0] = -params_.system_radius;
        pos[1] = -params_.system_radius;
        double d = 0.5 * (*spec)->GetSpecDiameter();
        double l = 0.25 * (*spec)->GetSpecLength();
        int num_x = (int)floor(2 * params_.system_radius / d);
        int num_y = (int)floor(2 * params_.system_radius / l);
        std::vector<std::pair<int, int>> grid_array;
        for (int i = 0; i < num_x; ++i) {
          for (int j = 0; j < num_y; ++j) {
            grid_array.push_back(std::make_pair(i, j));
          }
        }
        int *grid_index = new int[num_x * num_y];
        for (int i = 0; i < num_x * num_y; ++i) {
          grid_index[i] = i;
        }
        gsl_ran_shuffle(rng_.r, grid_index, num_x * num_y, sizeof(int));
        for (int i = 0; i < num_x * num_y; ++i) {
          pos[0] = grid_array[grid_index[i]].first * d;
          pos[1] = grid_array[grid_index[i]].second * l;
          (*spec)->AddMember();
          (*spec)->SetLastMemberPosition(pos);
          // First check that we are respecting boundary conditions
          std::vector<Object *> last_ixors;
          (*spec)->GetLastInteractors(&last_ixors);
          if (params_.boundary != 0 && !processing &&
              iengine_.CheckBoundaryConditions(last_ixors)) {
            (*spec)->PopMember();
            // We are not counting boundary condition failures in insertion
            // failures, since insertion failures are for packing issues
          }
          // Check if we have an overlap of objects
          else if (!force_overlap && !(*spec)->CanOverlap() && !processing &&
                   iengine_.CheckOverlap(last_ixors)) {
            (*spec)->PopMember();
            num_failures++;
          }
          // Otherwise update display of percentage of species inserted
          else {
            inserted++;
            iengine_.AddInteractors(last_ixors);
            int insert_percentage = (int)(100 * (float)inserted / (float)num);
            if (params_.print_complete) {
              if (insert_percentage % 10 == 0) {
                printf("  Inserting species: %d%% complete\n",
                       insert_percentage);
              }
            } else {
              printf("\r  Inserting species: %d%% complete", insert_percentage);
              fflush(stdout);
            }
          }
          /* Update the number of particles we need to insert, in case a species
             needs to have a certain packing fraction */
          num = (*spec)->GetNInsert();
          if (inserted == num)
            break;
        }
        delete[] grid_index;
      }
      putchar('\n');
      if (num != inserted) {
        printf("  Species insertion failure threshold of %d reached. "
               "Reattempting insertion.\n",
               params_.species_insertion_failure_threshold);
        (*spec)->PopAll();
        iengine_.Reset();
      }
      if (++num_attempts > params_.species_insertion_reattempt_threshold) {
        error_exit("Unable to insert species randomly within the reattempt "
                   "threshold of %d.\n",
                   params_.species_insertion_reattempt_threshold);
      }
    }
    if (!processing) {
      printf("\n");
      if ((*spec)->GetInsertionType().find("random") == std::string::npos) {
        (*spec)->ArrangeMembers();
        /*
         * This catch is breaking centered_oriented insertion. Skip it for now.
         *
         *if (!(*spec)->CanOverlap() &&
         *iengine_.CheckOverlap((*spec)->GetLastInteractors())) {
         *  error_exit("Species inserted with deterministic insertion type is
         *overlapping!");
         *}
         *
         */
      }
    }
  }
  /* Should do this all the time to force object counting */
  if (!processing) {
    iengine_.CheckUpdateObjects(); // Forces update as well
  }
}

/* Tear down data structures, e.g. cell lists, and close graphics window if
 * necessary. */
void Simulation::ClearSimulation() {
  output_mgr_.Close();
  ClearSpecies();
  iengine_.Clear();
#ifndef NOGRAPH
  if (params_.graph_flag) {
    graphics_.Clear();
  }
#endif
  std::cout << "  Simulation complete" << std::endl;
}

/* Tear down object type data structures */
void Simulation::ClearSpecies() {
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  species_factory_.clear();
}

/* Update the OpenGL graphics window */
void Simulation::Draw(bool single_frame) {
#ifndef NOGRAPH
  if (params_.graph_flag && i_step_ % params_.n_graph == 0) {
    /* Get updated object positions and orientations */
    GetGraphicsStructure();
    graphics_.Draw();
    if (params_.movie_flag) {
      /* Record bmp image of frame */
      grabber(graphics_.windx_, graphics_.windy_, params_.movie_directory,
              (single_frame ? 0 : frame_num_++));
    }
  }
#endif
}

/* Loop through objects in simulation to update their positions and
 * orientations in the graphics window */
void Simulation::GetGraphicsStructure() {
  graph_array_.clear();
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->Draw(&graph_array_);
  }
  /* Visualize interaction forces, crosslinks, etc */
  iengine_.DrawInteractions(&graph_array_);
}

/* Initialize output files */
void Simulation::InitOutputs() {
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), &i_step_,
                   run_name_);
  iengine_.InitOutputs();
  /* If analyzing run time, record cpu time here */
  if (params_.time_analysis) {
    cpu_init_time_ = cpu_time();
  }
}

/* Determine which output files we are reading */
void Simulation::InitInputs(run_options run_opts) {
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), &i_step_, run_name_,
                   true, run_opts.use_posits, run_opts.with_reloads,
                   run_opts.reduce_flag, run_opts.reduce_factor);
  iengine_.InitOutputs(true, run_opts.with_reloads, run_opts.reduce_flag);
}

/* Write object positions, etc if necessary */
void Simulation::WriteOutputs() {
  output_mgr_.WriteOutputs();
  /* Write interaction information/crosslink positions, etc */
  iengine_.WriteOutputs();
  /* If we are analyzing run time and this is the last step, record final time
   * here. */
  if (params_.time_analysis && i_step_ == params_.n_steps) {
    double cpu_final_time = cpu_time();
    double cpu_time = cpu_final_time - cpu_init_time_;
    std::cout << "CPU Time for Initialization: " << cpu_init_time_ << "\n";
    std::cout << "CPU Time: " << cpu_time << "\n";
    std::cout << "Sim Time: " << time_ << "\n";
    std::cout << "CPU Time/Sim Time: "
              << "\n"
              << cpu_time / time_ << std::endl;
  }
}

/* Run the steps we need for post-processing output files */
void Simulation::ProcessOutputs(system_parameters params,
                                run_options run_opts) {
  params_ = params;
  run_name_ = params.run_name;
  // Ensure that we are not trying to load any checkpoints when processing
  // outputs
  params_.load_checkpoint = 0;
  InitProcessing(run_opts);
  RunProcessing(run_opts);
  ClearSimulation();
}

// Initialize data structures for post-processing
void Simulation::InitProcessing(run_options run_opts) {
  space_.Init(&params_);
  InitObjects();
  InitSpecies();
  InsertSpecies(true, true);
  // if (run_opts.analysis_flag) {
  iengine_.Init(&params_, &species_, space_.GetStruct(), &i_step_, true);
  //}
  InitInputs(run_opts);
  if (run_opts.graphics_flag || run_opts.make_movie) {
    params_.graph_flag = 1;
    if (run_opts.use_posits && params_.n_graph < output_mgr_.GetNPosit()) {
      params_.n_graph = output_mgr_.GetNPosit();
    } else if (!run_opts.use_posits &&
               params_.n_graph < output_mgr_.GetNSpec()) {
      params_.n_graph = output_mgr_.GetNSpec();
    }
    if (run_opts.make_movie) {
      params_.movie_flag = 1;
    }
    InitGraphics();
  } else {
    params_.graph_flag = 0;
  }
  // if (run_opts.analysis_flag) {
  // for (auto it=species_.begin(); it!=species_.end(); ++it) {
  //(*it)->InitAnalysis();
  //}
  //}
}

/* Post-processing on simulation outputs for movie generation, analysis output
 * generation, etc. */
void Simulation::RunProcessing(run_options run_opts) {
  std::cout << "Processing outputs for: " << run_name_ << std::endl;
  bool local_order =
      (params_.local_order_analysis || params_.polar_order_analysis ||
       params_.overlap_analysis || params_.density_analysis);
  // Only step to n_steps-1 since we already read in one input at initialization
  int last_step =
      (run_opts.with_reloads ? params_.n_steps - 1 : 2 * params_.n_steps);
  bool run_analyses = run_opts.analysis_flag;
  for (i_step_ = 1; true; ++i_step_) {
    params_.i_step = i_step_;
    // i_step_ = params_.i_step_;
    time_ = (i_step_)*params_.delta;
    PrintComplete();
    output_mgr_.ReadInputs();
    iengine_.ReadInputs();
    if (early_exit) {
      early_exit = false;
      std::cout << "  Early exit triggered. Ending simulation.\n";
      if (run_analyses && i_step_ > params_.n_steps_equil) {
        for (auto it = species_.begin(); it != species_.end(); ++it) {
          (*it)->FinalizeAnalysis();
        }
      }
      return;
    }
    if (i_step_ <= params_.n_steps_equil) {
      Draw(run_opts.single_frame);
      continue;
    } else if (i_step_ == params_.n_steps_equil + 1 && run_analyses) {
      // InitAnalysis initalizes and runs the first batch of analyses
      for (auto it = species_.begin(); it != species_.end(); ++it) {
        (*it)->InitAnalysis();
      }
      // if (local_order && i_step_ % params_.local_order_n_analysis == 0) {
      // iengine_.StructureAnalysis();
      //}
      continue;
    }
    if (run_analyses) {
      bool struct_update = false;
      /* Check if we are running any species analysis to determine whether we
       * run structure analysis */
      for (auto it = species_.begin(); it != species_.end(); ++it) {
        if (((*it)->GetPositFlag() && i_step_ % (*it)->GetNPosit() == 0) ||
            ((*it)->GetSpecFlag() && i_step_ % (*it)->GetNSpec() == 0)) {
          struct_update = true;
        }
      }
      // Do structure analysis first
      if (struct_update && local_order &&
          i_step_ % params_.local_order_n_analysis == 0) {
        iengine_.StructureAnalysis();
      }
      // Now do species analyses
      for (auto it = species_.begin(); it != species_.end(); ++it) {
        if (((*it)->GetPositFlag() && i_step_ % (*it)->GetNPosit() == 0) ||
            ((*it)->GetSpecFlag() && i_step_ % (*it)->GetNSpec() == 0)) {
          (*it)->RunAnalysis();
        }
      }
    }
    Draw(run_opts.single_frame);
  }
  if (run_analyses) {
    for (auto it = species_.begin(); it != species_.end(); ++it) {
      (*it)->FinalizeAnalysis();
    }
  }
}
