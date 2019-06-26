#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_

class species_parameters;

class species_parameters {
  public:
    int num = 0;
    std::string insertion_type = "random";
    std::string insert_file = "none";
    int overlap = 0;
    std::string draw_type = "orientation";
    double color = 0;
    int posit_flag = 0;
    int spec_flag = 0;
    int checkpoint_flag = 0;
    int n_posit = 100;
    int n_spec = 100;
    int n_checkpoint = 10000;
};

class filament_parameters : public species_parameters {
  public:
    double diameter = 1;
    double length = -1;
    double persistence_length = 400;
    double max_length = 500;
    double min_length = 1;
    double max_bond_length = 5;
    double min_bond_length = 1.5;
    int spiral_flag = 0;
    double spiral_number_fail_condition = 0;
    double driving_factor = 0;
    double friction_ratio = 2;
    int dynamic_instability_flag = 0;
    int force_induced_catastrophe_flag = 0;
    int optical_trap_flag = 0;
    double optical_trap_spring = 20;
    int cilia_trap_flag = 0;
    double fic_factor = 0.828;
    double f_shrink_to_grow = 0.017;
    double f_shrink_to_pause = 0.0;
    double f_pause_to_grow = 0.0;
    double f_pause_to_shrink = 0.0;
    double f_grow_to_pause = 0.0;
    double f_grow_to_shrink = 0.00554;
    int metric_forces = 1;
    double v_poly = 0.44;
    double v_depoly = 0.793;
    int theta_analysis = 0;
    int lp_analysis = 0;
    int global_order_analysis = 0;
    double packing_fraction = -1;
    int shuffle = 1;
    double shuffle_factor = 10;
    double shuffle_frequency = 0.001;
    double perlen_ratio = -1;
    double perlen_ratio_target = -1;
    int n_bonds = -1;
    int driving_method = 0;
    int n_equil = 0;
    int orientation_corr_analysis = 0;
    int orientation_corr_n_steps = 1000;
    int crossing_analysis = 0;
    double intrinsic_curvature = 0;
    int flagella_flag = 0;
    double flagella_freq = 1;
    double flagella_period = 2;
    double flagella_amplitude = 1;
    int flocking_analysis = 0;
};

class passive_filament_parameters : public species_parameters {
  public:
    double diameter = 1;
    double length = -1;
    double persistence_length = 400;
    double max_length = 500;
    double min_length = 1;
    double max_bond_length = 5;
    double min_bond_length = 1.5;
    double driving_factor = 0;
    double friction_ratio = 2;
    int metric_forces = 1;
    int theta_analysis = 0;
    int lp_analysis = 0;
    double packing_fraction = -1;
    double perlen_ratio = -1;
    int n_bonds = -1;
};

class hard_rod_parameters : public species_parameters {
  public:
    double diameter = 1;
    double length = 40;
    double min_length = 3;
    double max_length = 300;
    double max_bond_length = 5;
    double driving_factor = 0;
};

class br_bead_parameters : public species_parameters {
  public:
    double diameter = 1;
    double driving_factor = 0;
    double packing_fraction = -1;
};

class md_bead_parameters : public species_parameters {
  public:
    double diameter = 1;
    double mass = 1;
    double driving_factor = 0;
};

class centrosome_parameters : public species_parameters {
  public:
    double diameter = 10;
    int n_filaments_min = 0;
    int n_filaments_max = 10;
    int fixed_spacing = 0;
    int alignment_potential = 0;
    double k_spring = 1000;
    double k_align = 0;
    double spring_length = 0;
};

class motor_parameters : public species_parameters {
  public:
    double diameter = 3;
    int walker = 0;
    int step_direction = 0;
    double k_off = 2;
    double k_on = 10;
    double concentration = 0;
    double velocity = 1;
    int diffusion_flag = 1;
    double k_spring = 1000;
    double f_spring_max = 100;
};

class bead_spring_parameters : public species_parameters {
  public:
    double diameter = 1;
    double length = 40;
    double persistence_length = 4000;
    double max_bond_length = 1;
    double bond_rest_length = 0.8;
    double bond_spring = 100;
    double driving_factor = 0;
    int lp_analysis = 0;
    int theta_analysis = 0;
    double packing_fraction = -1;
};

class spherocylinder_parameters : public species_parameters {
  public:
    double diameter = 1;
    double length = 5;
    int diffusion_analysis = 0;
    int n_diffusion_samples = 1;
    int midstep = 0;
};

class spindle_parameters : public species_parameters {
  public:
    double diameter = 10;
    double length = 20;
    int n_filaments_bud = 1;
    int n_filaments_mother = 0;
    int alignment_potential = 0;
    double k_spring = 1000;
    double k_align = 0;
    double spring_length = 0;
    double spb_diameter = 5;
};

class crosslink_parameters : public species_parameters {
  public:
    double concentration = 0;
    double diameter = 1;
    int walker = 0;
    double velocity = 1;
    int diffusion_flag = 0;
    double k_on = 10;
    double k_off = 2;
    double k_on_sd = 10;
    double force_dep_factor = 1;
    double polar_affinity = 0;
    double k_spring = 10;
    double f_spring_max = 1000;
    double f_stall = 100;
    int force_dep_vel_flag = 1;
    double k_align = 0;
    double rest_length = 0;
    int step_direction = 0;
    std::string tether_draw_type = "orientation";
    double tether_diameter = 0.5;
    double tether_color = 3.1416;
    int end_pausing = 0;
    double r_capture = 5;
};

class system_parameters {
  public:
    std::string default_param_file = "src/config_params.yaml";
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
    double cell_length = 10;
    int n_update_cells = 0;
    int graph_flag = 0;
    int n_graph = 1000;
    double graph_diameter = 0;
    int graph_background = 1;
    int draw_boundary = 1;
    int load_checkpoint = 0;
    std::string checkpoint_run_name = "sc";
    int n_load = 0;
    int print_complete = 0;
    std::string insertion_type = "species";
    int movie_flag = 0;
    std::string movie_directory = "frames";
    int time_analysis = 0;
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
    int constant_pressure = 0;
    int constant_volume = 0;
    double target_pressure = 0;
    double target_radius = 100;
    int pressure_time = 100;
    double compressibility = 1;
    int stoch_flag = 1;
    int thermo_flag = 0;
    int n_thermo = 1000;
    int interaction_flag = 1;
    int species_insertion_failure_threshold = 10000;
    int species_insertion_reattempt_threshold = 10;
    int uniform_crystal = 0;
    int n_steps_equil = 0;
    int n_steps_target = 100000;
    int static_particle_number = 0;
    int checkpoint_from_spec = 0;
    std::string potential = "wca";
    double soft_potential_mag = 10;
    double soft_potential_mag_target = -1;
    int like_like_interactions = 1;
    int auto_graph = 0;
    double local_order_width = 50;
    double local_order_bin_width = 0.5;
    int local_order_average = 1;
    int local_order_analysis = 0;
    int local_order_n_analysis = 100;
    int density_analysis = 0;
    double density_bin_width = 0.1;
    int density_com_only = 0;
    int polar_order_analysis = 0;
    int polar_order_n_bins = 100;
    double polar_order_contact_cutoff = 3;
    int polar_order_color = 0;
    int overlap_analysis = 0;
    int highlight_overlaps = 0;
    int reduced = 0;
    int reload_reduce_switch = 0;
    double flock_polar_min = 0.5;
    double flock_contact_min = 1.5;
    int highlight_flock = 0;
    double flock_color_ext = 1.57;
    double flock_color_int = 4.71;
    int in_out_flag = 0;
    species_parameters species;
    filament_parameters filament;
    passive_filament_parameters passive_filament;
    hard_rod_parameters hard_rod;
    br_bead_parameters br_bead;
    md_bead_parameters md_bead;
    centrosome_parameters centrosome;
    motor_parameters motor;
    bead_spring_parameters bead_spring;
    spherocylinder_parameters spherocylinder;
    spindle_parameters spindle;
    crosslink_parameters crosslink;
};

#endif // _SIMCORE_PARAMETERS_H_