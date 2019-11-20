#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_

#include "definitions.hpp"

#include <string>

template <unsigned char S> struct species_parameters {
  std::string name = "species";
  int num = 0;
  double diameter = 1;
  double length = 0;
  std::string insertion_type = "random";
  std::string insert_file = "none";
  bool overlap = false;
  std::string draw_type = "orientation";
  double color = 0;
  bool posit_flag = false;
  bool spec_flag = false;
  int n_posit = 100;
  int n_spec = 100;
  virtual ~species_parameters() {}
};

typedef species_parameters<species_id::none> species_base_parameters;

template <>
struct species_parameters<species_id::rigid_filament>
    : public species_base_parameters {
  double max_length = 500;
  double min_length = 5;
  bool stationary_flag = false;
  double packing_fraction = -1;
  int n_equil = 0;
};
typedef species_parameters<species_id::rigid_filament> rigid_filament_parameters;

template <>
struct species_parameters<species_id::filament>
    : public species_base_parameters {
  double persistence_length = 400;
  double max_length = 500;
  double min_length = 5;
  double min_bond_length = 1.5;
  bool spiral_flag = false;
  double spiral_number_fail_condition = 0;
  double driving_factor = 0;
  double friction_ratio = 2;
  bool dynamic_instability_flag = false;
  bool force_induced_catastrophe_flag = false;
  bool optical_trap_flag = false;
  double optical_trap_spring = 20;
  bool optical_trap_fixed = false;
  bool cilia_trap_flag = false;
  double fic_factor = 0.828;
  double f_shrink_to_grow = 0.017;
  double f_shrink_to_pause = 0.0;
  double f_pause_to_grow = 0.0;
  double f_pause_to_shrink = 0.0;
  double f_grow_to_pause = 0.0;
  double f_grow_to_shrink = 0.00554;
  double v_poly = 0.44;
  double v_depoly = 0.793;
  bool theta_analysis = false;
  bool lp_analysis = false;
  bool global_order_analysis = false;
  double packing_fraction = -1;
  double perlen_ratio = -1;
  bool drive_from_bond_center = true;
  int n_equil = 0;
  bool orientation_corr_analysis = false;
  int orientation_corr_n_steps = 1000;
  bool crossing_analysis = false;
  double intrinsic_curvature = 0;
  bool flagella_flag = false;
  double flagella_freq = 1;
  double flagella_period = 2;
  double flagella_amplitude = 1;
  bool flocking_analysis = false;
  bool polydispersity_flag = false;
};
typedef species_parameters<species_id::filament> filament_parameters;

template <>
struct species_parameters<species_id::br_bead>
    : public species_base_parameters {
  double driving_factor = 0;
  double packing_fraction = -1;
};
typedef species_parameters<species_id::br_bead> br_bead_parameters;

template <>
struct species_parameters<species_id::spherocylinder>
    : public species_base_parameters {
  bool diffusion_analysis = false;
  int n_diffusion_samples = 1;
  bool midstep = false;
};
typedef species_parameters<species_id::spherocylinder> spherocylinder_parameters;

template <>
struct species_parameters<species_id::spindle>
    : public species_base_parameters {
  int n_filaments_bud = 1;
  int n_filaments_mother = 0;
  bool alignment_potential = false;
  double k_spring = 1000;
  double k_align = 0;
  double spring_length = 0;
  double spb_diameter = 5;
};
typedef species_parameters<species_id::spindle> spindle_parameters;

template <>
struct species_parameters<species_id::crosslink>
    : public species_base_parameters {
  double concentration = 0;
  bool walker_flag = false;
  bool static_flag = false;
  bool diffusion_flag = false;
  double velocity = 1;
  double k_on = 10;
  double k_off = 2;
  double k_on_d = 10;
  double k_off_d = 2;
  double force_dep_factor = 1;
  double polar_affinity = 0;
  double k_spring = 10;
  double f_stall = 100;
  bool force_dep_vel_flag = true;
  double k_align = 0;
  double rest_length = 0;
  int step_direction = 0;
  std::string tether_draw_type = "orientation";
  double tether_diameter = 0.5;
  double tether_color = 3.1416;
  bool end_pausing = false;
  double r_capture = 5;
};
typedef species_parameters<species_id::crosslink> crosslink_parameters;

struct system_parameters {
  long seed = 7859459105545;
  int n_runs = 1;
  int n_random = 1;
  std::string run_name = "sc";
  int n_dim = 3;
  int n_periodic = 0;
  int boundary = 0;
  double system_radius = 100;
  int n_steps = 1000000;
  int i_step = 0;
  double delta = 0.001;
  bool graph_flag = false;
  int n_graph = 1000;
  double graph_diameter = 0;
  bool invert_background = false;
  bool draw_boundary = true;
  bool load_checkpoint = false;
  std::string checkpoint_run_name = "sc";
  int n_load = 0;
  bool movie_flag = false;
  std::string movie_directory = "frames";
  bool time_analysis = false;
  double bud_height = 680;
  double bud_radius = 300;
  double lj_epsilon = 1;
  double wca_eps = 1;
  double wca_sig = 1;
  double ss_a = 1;
  double ss_rs = 1.5;
  double ss_eps = 1;
  double f_cutoff = 2000;
  bool constant_pressure = false;
  bool constant_volume = false;
  double target_pressure = 0;
  double target_radius = 100;
  int pressure_time = 100;
  double compressibility = 1;
  bool stoch_flag = true;
  bool thermo_flag = false;
  int n_thermo = 1000;
  bool interaction_flag = true;
  int species_insertion_failure_threshold = 10000;
  int species_insertion_reattempt_threshold = 10;
  bool uniform_crystal = false;
  int n_steps_equil = 0;
  int n_steps_target = 100000;
  bool static_particle_number = false;
  bool checkpoint_from_spec = false;
  std::string potential = "wca";
  double soft_potential_mag = 10;
  double soft_potential_mag_target = -1;
  bool like_like_interactions = true;
  bool auto_graph = false;
  bool local_order_analysis = false;
  double local_order_width = 50;
  double local_order_bin_width = 0.5;
  int local_order_n_analysis = 100;
  int density_analysis = 0;
  double density_bin_width = 0.1;
  bool density_com_only = false;
  bool polar_order_analysis = false;
  int polar_order_n_bins = 100;
  double polar_order_contact_cutoff = 3;
  bool overlap_analysis = false;
  bool highlight_overlaps = false;
  bool reduced = false;
  bool reload_reduce_switch = false;
  double flock_polar_min = 0.5;
  double flock_contact_min = 1.5;
  bool highlight_flock = false;
  double flock_color_ext = 1.57;
  double flock_color_int = 4.71;
  bool in_out_flag = false;
  bool checkpoint_flag = false;
  int n_checkpoint = 10000;
  bool no_midstep = false;
};

#endif // _SIMCORE_PARAMETERS_H_