#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_

struct filament_parameters {
  int overlap = 0;
  std::string insertion_type = "random";
  double max_child_length = 3;
  int num = 0;
  double color = 0;
  std::string draw_type = "orientation";
  double diameter = 1;
  double length = 40;
  double persistence_length = 400;
  double max_length = 300;
  double min_length = 5;
  int spiral_flag = 0;
  double driving_factor = 0;
  double friction_ratio = 2;
  int dynamic_instability_flag = 0;
  int force_induced_catastrophe_flag = 0;
  double f_shrink_to_grow = 0.017;
  double f_shrink_to_pause = 0.0;
  double f_pause_to_grow = 0.0;
  double f_pause_to_shrink = 0.0;
  double f_grow_to_pause = 0.0;
  double f_grow_to_shrink = 0.00554;
  int metric_forces = 1;
  double v_poly = 0.44;
  double v_depoly = 0.793;
  int posit_flag = 0;
  int n_posit = 100;
  std::string posit_file = "filament.posit";
};

struct hard_rod_parameters {
  std::string insertion_type = "random";
  int overlap = 0;
  int num = 0;
  double color = 0;
  std::string draw_type = "orientation";
  double diameter = 1;
  double length = 40;
  double min_length = 3;
  double max_length = 300;
  double max_child_length = 5;
  double driving_factor = 0;
  int posit_flag = 0;
  int n_posit = 100;
  std::string posit_file = "hard_rod.posit";
};

struct br_bead_parameters {
  std::string insertion_type = "random";
  int overlap = 0;
  int num = 0;
  double color = 0;
  std::string draw_type = "flat";
  double diameter = 1;
  int posit_flag = 0;
  double driving_factor = 0;
  int n_posit = 100;
  std::string posit_file = "br_bead.posit";
};

struct md_bead_parameters {
  std::string insertion_type = "random";
  int overlap = 0;
  int num = 0;
  double color = 0;
  std::string draw_type = "flat";
  double diameter = 1;
  double mass = 1;
  int posit_flag = 0;
  int n_posit = 100;
  std::string posit_file = "md_bead.posit";
};

struct system_parameters {
  long seed = 7859459105545;
  int n_runs = 1;
  int n_random = 1;
  std::string run_name = "sc";
  int n_dim = 3;
  int n_periodic = 0;
  int boundary_type = 0;
  double system_radius = 100;
  int n_steps = 1000000;
  double delta = 0.001;
  double cell_length = 10;
  int n_update_cells = 10;
  int graph_flag = 0;
  int n_graph = 1000;
  double graph_diameter = 0;
  int graph_background = 1;
  int movie_flag = 0;
  std::string movie_directory = "frames";
  int time_flag = 0;
  double bud_height = 680;
  double bud_radius = 300;
  double lj_epsilon = 1;
  double wca_eps = 1;
  double wca_sig = 1;
  double ss_a = 1;
  double ss_rs = 1.5;
  double ss_eps = 1;
  double f_cutoff = 100;
  int max_overlap = 100000;
  int virial_time_avg = 100;
  int constant_pressure = 0;
  int constant_volume = 0;
  double target_pressure = 0;
  double target_radius = 100;
  int pressure_time = 100;
  double compressibility = 1;
  int posit_flag = 0;
  int n_posit = 100;
  filament_parameters filament;
  hard_rod_parameters hard_rod;
  br_bead_parameters br_bead;
  md_bead_parameters md_bead;
};

#endif // _SIMCORE_PARAMETERS_H_