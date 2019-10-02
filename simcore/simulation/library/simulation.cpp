
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
  Logger::Info("Running simulation on %d threads", omp_get_max_threads());
#else
  Logger::Info("Running simulation");
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
      Logger::Info("Early exit triggered. Ending simulation.");
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
  long iteration = i_step_ * 10;
  long steps = params_.n_steps;
  if (iteration % steps == 0) {
    Logger::Info("%d%% complete", 10 * iteration / steps);
  }
  Logger::Trace("*****Step %d*****", i_step_);
}

/* Update the positions of all objects in the system using numerical
 * integration for one time step defined by the delta parameter. */
void Simulation::Integrate() {
  Logger::Debug("Updating object positions");
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->UpdatePositions();
  }
}

/* Calculate interaction forces between all objects if necessary. */
void Simulation::Interact() {
  Logger::Debug("Calculating object interactions");
  iengine_.Interact();
}

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
  Logger::Debug("Calculating system thermodynamics");
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
  Logger::Info("Initializing simulation");
#ifdef TRACE
  if (params_.n_steps > 100) {
    Logger::Warning("Simulation run in trace mode with a large number of "
                    "n_steps (%d).",
                    params_.n_steps);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif
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
#ifdef TRACE
      if (spec->GetNInsert() > 20) {
        Logger::Warning("Simulation run in trace mode with a large number of "
                        "objects in species %s (%d).",
                        spec->GetSID()._to_string(), spec->GetNInsert());
        fprintf(stderr, "Continue anyway? (y/N) ");
        char c;
        if (std::cin.peek() != 'y') {
          Logger::Error("Terminating simulation by user request");
        } else if (!(std::cin >> c)) {
          Logger::Error("Invalid input");
        } else {
          fprintf(stderr, "Resuming simulation\n");
          std::cin.ignore();
        }
      }
#endif
      species_.push_back(spec);
      species_.back()->Reserve();
    }
  }
}

/* Initialize object positions and orientations.*/
void Simulation::InsertSpecies(bool force_overlap, bool processing) {
  Logger::Info("Inserting species");
  for (auto spec = species_.begin(); spec != species_.end(); ++spec) {
    // Check for random insertion
    if (processing ||
        (*spec)->GetInsertionType().find("random") == std::string::npos) {
      /* Insertion is non-random: don't check for overlaps */
      force_overlap = true;
    }
    int num = (*spec)->GetNInsert();
    bool not_done = true;
    int inserted = 0;
    int num_attempts = 0;
    while (num != inserted) {
      inserted = 0;
      int num_failures = 0;
      while (num != inserted) {
        (*spec)->AddMember();
        /* Update the number of particles we need to insert, in case a species
           needs to have a certain packing fraction */
        num = (*spec)->GetNInsert();
        // First check that we are respecting boundary conditions
        std::vector<Object *> last_ixors;
        (*spec)->GetLastInteractors(&last_ixors);
        if (params_.boundary != 0 && !processing &&
            iengine_.CheckBoundaryConditions(last_ixors)) {
          (*spec)->PopMember();
          /* We are not counting boundary condition failures in insertion
           failures, since insertion failures are for packing issues */
        }
        // Check if we have an overlap of objects
        else if (!force_overlap && !(*spec)->CanOverlap() && !processing &&
                 iengine_.CheckOverlap(last_ixors)) {
          (*spec)->PopMember();
          num_failures++;
        }
        /* Otherwise update interaction engine to include new interactors and
           update display of percentage of species inserted */
        else {
          inserted++;
          if (!force_overlap && !processing) {
            iengine_.AddInteractors(last_ixors);
          }
        }
        if (num_failures > params_.species_insertion_failure_threshold) {
          Logger::Warning("Too many insertion failures have occurred. Managed "
                          "to insert %2.2f% of objects",
                          100.0 * inserted / (float)num);
          break;
        }
      }
      if (num != inserted) {
        // Attempt a lattice-based insertion strategy (only 2d for now)
        if (params_.n_dim == 3)
          continue;
        Logger::Warning("Attempting lattice-based insertion strategy");
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
          /* Update the number of particles we need to insert, in case a
             species needs to have a certain packing fraction */
          num = (*spec)->GetNInsert();
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
          }
          if (inserted == num)
            break;
        }
        delete[] grid_index;
      }
      if (num != inserted) {
        Logger::Warning("Species insertion failure threshold of %d reached. "
                        "Reattempting insertion.\n",
                        params_.species_insertion_failure_threshold);
        (*spec)->PopAll();
        iengine_.Reset();
      }
      if (++num_attempts > params_.species_insertion_reattempt_threshold) {
        Logger::Error("Unable to insert species randomly within the reattempt "
                      "threshold of %d.\n",
                      params_.species_insertion_reattempt_threshold);
      }
    }
    if (!processing) {
      if ((*spec)->GetInsertionType().find("random") == std::string::npos) {
        (*spec)->ArrangeMembers();
      }
    }
    if (params_.load_checkpoint) {
      (*spec)->LoadFromCheckpoints(run_name_, params_.checkpoint_run_name);
    }
  }
  /* Should do this all the time to force object counting */
  iengine_.ResetCellList(); // Forces rebuild cell list without redundancy
  iengine_.Reset();
  //if (!processing) {
  iengine_.CheckUpdateObjects(); // Forces update as well
  //}
}

/* Tear down data structures, e.g. cell lists, and close graphics window if
 * necessary. */
void Simulation::ClearSimulation() {
  Logger::Debug("Clearing simulation resources");
  output_mgr_.Close();
  ClearSpecies();
  iengine_.Clear();
#ifndef NOGRAPH
  if (params_.graph_flag) {
    graphics_.Clear();
  }
#endif
  Logger::Info("Simulation complete");
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
  Logger::Trace("Drawing graphable objects");
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
  Logger::Trace("Retrieving graphics structures from objects");
  graph_array_.clear();
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->Draw(&graph_array_);
  }
  /* Visualize interaction forces, crosslinks, etc */
  iengine_.DrawInteractions(&graph_array_);
}

/* Initialize output files */
void Simulation::InitOutputs() {
  Logger::Debug("Initializing output files");
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), &i_step_,
                   run_name_);
  //if (!params_.load_checkpoint)
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
  Logger::Debug("Writing outputs");
  output_mgr_.WriteOutputs();
  /* Write interaction information/crosslink positions, etc */
  iengine_.WriteOutputs();
  /* If we are analyzing run time and this is the last step, record final time
   * here. */
  if (params_.time_analysis && i_step_ == params_.n_steps) {
    double cpu_final_time = cpu_time();
    double cpu_time = cpu_final_time - cpu_init_time_;
    Logger::Info("CPU Time for Initialization: %2.6f", cpu_init_time_);
    Logger::Info("CPU Time: %2.6f", cpu_time);
    Logger::Info("Sim Time: %2.6f", time_);
    Logger::Info("CPU Time/Sim Time: %2.6f\n", cpu_time / time_);
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
  Logger::Info("Initializing datastructures for post-processing outputs");
  space_.Init(&params_);
  InitObjects();
  InitSpecies();
  iengine_.Init(&params_, &species_, space_.GetStruct(), &i_step_, true);
  InsertSpecies(true, true);
  // if (run_opts.analysis_flag) {
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
  Logger::Info("Processing outputs for %s", run_name_.c_str());
  bool local_order =
      (params_.local_order_analysis || params_.polar_order_analysis ||
       params_.overlap_analysis || params_.density_analysis);
  // Only step to n_steps-1 since we already read in one input at
  // initialization
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
      Logger::Info("Early exit triggered. Ending simulation.");
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
