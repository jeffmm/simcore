#include "simcore/filament_species.hpp"

FilamentSpecies::FilamentSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::filament);
}
void FilamentSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  fill_volume_ = 0;
  packing_fraction_ = sparams_.packing_fraction;
#ifdef TRACE
  if (packing_fraction_ > 0) {
    Logger::Warning("Simulation run in trace mode with a potentially large "
                    "number of objects in species %s (packing fraction ="
                    " %2.4f)",
                    GetSID()._to_string(), packing_fraction_);
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

  double min_length = 2 * sparams_.min_bond_length;
  if (!sparams_.polydispersity_flag && sparams_.length < min_length) {
    Logger::Warning("Filament length %2.2f is less than minimum filament length"
                    " %2.2f for minimum bond length %2.2f. Setting length to "
                    "minimum length.",
                    sparams_.length, min_length, sparams_.min_bond_length);
    sparams_.length = min_length;
  }
  if (sparams_.perlen_ratio > 0) {
    if (sparams_.polydispersity_flag) {
      Logger::Warning("Persistence length ratio parameter is set with "
                      "polydispersity. Ignoring perlen_ratio and using "
                      "persistence_length parameter instead.");
    } else {
      sparams_.persistence_length = sparams_.length * sparams_.perlen_ratio;
    }
  }
  if (sparams_.flexure_number >= 0) {
    if (sparams_.polydispersity_flag) {
      Logger::Warning("Flexure number parameter is set with "
                      "polydispersity. Ignoring flexure number and using "
                      "persistence_length and driving_factor parameters "
                      "instead.");
    } else {
      sparams_.peclet_number = sparams_.flexure_number *
                               sparams_.persistence_length / sparams_.length;
    }
  }
  if (sparams_.peclet_number >= 0) {
    if (sparams_.dynamic_instability_flag) {
      Logger::Error(
          "Dynamic instability and Peclet/Flexure numbers are not compatible");
    }
    if (sparams_.polydispersity_flag) {
      Logger::Warning("Peclet number parameter is set with "
                      "polydispersity. Ignoring Peclet number and using "
                      "driving_factor parameter instead.");
    } else {
      sparams_.driving_factor = sparams_.peclet_number / SQR(sparams_.length);
    }
  }
  if (sparams_.spiral_init_flag && sparams_.spiral_number_fail_condition <= 0) {
    Logger::Warning("Spiral simulation will not end on spiral failure due to"
                    " nonpositive spiral_number_fail_condition");
  }

  if (sparams_.error_analysis) {
    if (sparams_.reference_frame_flag) {
      Logger::Warning(
          "Filament errors will be inflated due to rotation of site "
          "positions into filament reference frame. Errors should be viewed "
          "only"
          " in this context.");
    }
    InitErrorAnalysis();
  }
}

void FilamentSpecies::CleanUp() {
  Species::CleanUp();
  if (error_file_.is_open()) {
    error_file_.close();
  }
}

void FilamentSpecies::UpdatePositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  filament_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  filament_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i)
      for (auto it = chunks[i].first; it != chunks[i].second; ++it)
        it->UpdatePosition(midstep_);
  }
#else
  for (filament_iterator it = members_.begin(); it != members_.end(); ++it)
    it->UpdatePosition(midstep_);
#endif

  midstep_ = !midstep_;
  if (sparams_.error_analysis) {
    RunErrorAnalysis();
  }
  if (sparams_.spiral_init_flag && sparams_.spiral_number_fail_condition >= 0) {
    int n_failed_spirals = 0;
    for (auto it = members_.begin(); it != members_.end(); ++it) {
      if (ABS(it->GetSpiralNumber()) < sparams_.spiral_number_fail_condition) {
        // Failed spiral
        n_failed_spirals++;
      }
    }
    // If all spirals have failed, end simulation early
    if (n_failed_spirals == n_members_) {
      early_exit = true;
    }
  }
}

void FilamentSpecies::InitErrorAnalysis() {
  std::string fname =
      params_->run_name + "_filament_" + sparams_.name + ".error.analysis";
  error_file_.open(fname, std::ios::out);
  if (!error_file_.is_open()) {
    Logger::Error("Filament error analysis file %s failed to open!",
                  fname.c_str());
  }
  error_file_ << "n_dim delta length persistence_length min_bond_length "
                 "driving_factor\n";
  error_file_ << params_->n_dim << " " << params_->delta << " "
              << sparams_.length << " " << sparams_.persistence_length << " "
              << sparams_.min_bond_length << " " << sparams_.driving_factor
              << "\n";
  error_file_ << "steps_to_renorm\n";
}

void FilamentSpecies::RunErrorAnalysis() {
  if (!error_file_.is_open()) {
    Logger::Error("Error analysis file failed to open in FilamentSpecies");
  }
  std::vector<int> error_rates;
  for (filament_iterator it = members_.begin(); it != members_.end(); ++it)
    it->GetErrorRates(error_rates);
  for (auto it = error_rates.begin(); it != error_rates.end(); ++it) {
    error_file_ << *it << "\n";
  }
}

void FilamentSpecies::Reserve() {
  int max_insert = GetNInsert();
  if (packing_fraction_ > 0) {
    double min_length = 2 * sparams_.min_bond_length;
    double diameter = sparams_.diameter;
    double min_vol = 0;
    if (params_->n_dim == 2) {
      min_vol = diameter * min_length + 0.25 * M_PI * diameter * diameter;
    }
    if (params_->n_dim == 3) {
      min_vol = 0.25 * M_PI * diameter * diameter * min_length +
                1.0 / 6.0 * M_PI * diameter * diameter * diameter;
    }
    max_insert = (int)ceil(packing_fraction_ * space_->volume / min_vol);
  }
  members_.reserve(max_insert);
  Logger::Debug("Reserving memory for %d members in FilamentSpecies",
                max_insert);
}

void FilamentSpecies::AddMember() {
  Species::AddMember();
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ += vol;
    /* if we are still short on volume for the target packing fraction, then
       request more members */
    if (fill_volume_ < packing_fraction_ * space_->volume &&
        members_.size() == sparams_.num) {
      sparams_.num++;
    }
  }
}

void FilamentSpecies::PopMember() {
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ -= vol;
  }
  Species::PopMember();
}

void FilamentSpecies::ResetPreviousPositions() {
  Species::ResetPreviousPositions();
  midstep_ = true;
}

const double FilamentSpecies::GetSpecLength() const {
  if (sparams_.dynamic_instability_flag) {
    return 2 * sparams_.min_bond_length;
  } else {
    return 1.5 * sparams_.min_bond_length;
  }
}
void FilamentSpecies::LoadAnalysis() {
  if (sparams_.curvature_cluster_analysis) {
    FilamentAnalysis *ccluster = new CurvatureClusterAnalysis;
    analysis_.push_back(ccluster);
  }
  if (sparams_.spiral_analysis) {
    FilamentAnalysis *spiral = new SpiralAnalysis;
    analysis_.push_back(spiral);
  }
  if (sparams_.theta_analysis) {
    FilamentAnalysis *theta = new AngleDistributionAnalysis;
    analysis_.push_back(theta);
  }
  /* PolarOrderAnalysis before FlockingAnalysis, so we can calculate polar order
     only once if both analyses are being done  */
  if (sparams_.polar_order_analysis) {
    FilamentAnalysis *polar_order = new PolarOrderAnalysis;
    analysis_.push_back(polar_order);
  }
  if (sparams_.flocking_analysis) {
    FilamentAnalysis *flock = new FlockingAnalysis;
    analysis_.push_back(flock);
  }
  if (sparams_.orientation_corr_analysis) {
    FilamentAnalysis *ocorr = new OrientationCorrelationAnalysis;
    analysis_.push_back(ocorr);
  }
  if (sparams_.global_order_analysis) {
    FilamentAnalysis *global_order = new GlobalOrderAnalysis;
    analysis_.push_back(global_order);
  }
  if (sparams_.number_fluctuation_analysis) {
    FilamentAnalysis *gnf = new GiantNumberFluctuationAnalysis;
    analysis_.push_back(gnf);
  }
  if (sparams_.lp_analysis) {
    FilamentAnalysis *mse2e = new EndToEndFluctuationAnalysis;
    analysis_.push_back(mse2e);
  }
  if (sparams_.in_out_analysis) {
    FilamentAnalysis *in_out = new IncomingOutgoingAngleAnalysis;
    analysis_.push_back(in_out);
  }
  if (sparams_.crossing_analysis) {
    FilamentAnalysis *crossing = new BarrierCrossingAnalysis;
    analysis_.push_back(crossing);
  }
}
