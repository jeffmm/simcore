
#include "simulation.h"

Simulation::Simulation() {}

Simulation::Simulation(std::string param_file, std::string run_name, int n_runs) {
  params.init();
  ParseParams(param_file);
  rng.init(params.seed);
  for (i_run=0; i_run<n_runs; ++i_run) {
    InitSimulation();
    RunSimulation();
    ClearSimulation();
  }
}

void Simulation::RunSimulation() {
  for (i_step=0; i_step<params.n_steps; ++i_step) {
    time = i_step * params.delta; 
    Integrate();
    Draw();
    WriteOutputs();
  }
}

void Simulation::Integrate() {
  for (auto sys=systems_.begin(); sys!=systems_.end(); ++sys)
    (*sys)->Integrate();
}

void Simulation::InitSimulation() {
  InitSpace();
  InitSystems();
  //InitIntegrator();
  if (params.graph_flag)
    InitGraphics();
}

void Simulation::InitSystems() {
  ObjectSystemBase *sys_ptr = new ParticleMDSystem(&params, &space, gsl_rng_get(rng.r));
  systems_.push_back(sys_ptr);
  delete sys_ptr;
  for (auto sys=systems_.begin(); sys!=systems_.end(); ++sys)
    (*sys)->InitElements();
}

void Simulation::ClearSimulation() {
  ClearSpace();
  if (params.graph_flag)
    ClearGraphics();
}

void Simulation::ClearSpace() {
  space.Clear();
}

void Simulation::InitSpace() {
  space.Init(&params, gsl_rng_get(rng.r));
}

//void Simulation::InitIntegrator() {
  //integrator.Init(&params, gsl_rng_get(rng.r));
//}

//void Simulation::ClearIntegrator() {
  //integrator.Clear();
//}

void Simulation::ParseParams(std::string param_file) {
 
  std::cout << "Reading parameters from " << param_file << std::endl;
  YAML::Node node = YAML::LoadFile(param_file);

  for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
    std::string param_name = it->first.as<std::string>();
    std::string param_value = it->second.as<std::string>();
  
    // check for all parameters
#include "parse_params_body.cpp"
  }
  std::cout << std::endl;
}

void Simulation::InitGraphics() {
  int n_dim = params.n_dim;
  int max_objects = 10000; // FIXME should compute this from parameters
  g_struct.h = space.GetUnitCell();
  g_struct.r = new double*[max_objects];
  for(int i=0; i<max_objects; ++i) 
    g_struct.r[i] = new double[n_dim]; 
  g_struct.u = new double*[max_objects];
  for(int i=0; i<max_objects; ++i) 
    g_struct.u[i] = new double[n_dim]; 
  g_struct.l = new double[max_objects];
  g_struct.diam = new double[max_objects];
  g_struct.m_rad = params.mother_radius;
  g_struct.d_rad = params.daughter_radius;
  g_struct.m_d_dist = params.mother_daughter_dist;
  
  GetGraphicsStructure();
  double background_color = (params.graph_background == 0 ? 0.1 : 1);
  graphics.Init(params.n_dim, g_struct.h, background_color);
  if (params.boundary_type == 0)
    graphics.SetBoundaryType("sphere");
  else if (params.boundary_type == 1)
    graphics.SetBoundaryType("cube");
  else if (params.boundary_type == 2)
    graphics.SetBoundaryType("snowman");
  graphics.DrawLoop(g_struct);
}

void Simulation::ClearGraphics() {
  int max_objects = 10000; //FIXME
  for(int i=0; i<max_objects; ++i) {
    delete[] g_struct.r;
    delete[] g_struct.u;
  }
  delete[] g_struct.l;
  delete[] g_struct.diam;
  graphics.Clear();
}

void Simulation::Draw() {
  GetGraphicsStructure();
  graphics.Draw(g_struct);
  // Record bmp image of frame 
  if (params.grab_flag) {
    grabber(graphics.windx_, graphics.windy_,
            params.grab_file, (int) i_step/params.n_graph);
  }
}

void Simulation::GetGraphicsStructure() {

  int tot_spheros = 0; 
  for (auto sys=systems_.begin(); sys!=systems_.end(); ++sys) {
    tot_spheros += (*sys)->GetNElements();
    (*sys)->GetGraphics(&g_struct);
  }
  g_struct.n_spheros = tot_spheros;
  if (space.GetType() == SNOWMAN) {
    double z_correct = 0.5 * (params.mother_daughter_dist + params.mother_radius + params.daughter_radius) - params.mother_radius;
    for (int i=0; i<tot_spheros; ++i)
      g_struct.r[i][params.n_dim-1] -= z_correct;
  }
}

void Simulation::InitOutputs() {
  if (params.time_analysis_flag) {
    cpu_init_time = cpu();
  }
}

void Simulation::WriteOutputs() {
  if (i_step==0) {
    InitOutputs();
    return;
  }
  if (params.time_analysis_flag && i_step == params.n_steps-1) {
    double cpu_time = cpu() - cpu_init_time;
    std::cout << "CPU Time for Initialization: " <<  cpu_init_time << std::endl;
    std::cout << "CPU Time: " << cpu_time << std::endl;
    std::cout << "Sim Time: " << time << std::endl;
    std::cout << "CPU Time/Sim Time: " << std::endl << cpu_time/time << std::endl;
  }
}

