#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_

class species_parameters;

class species_parameters {
  public:
    int num = 0;
    std::string insertion_type = "random";
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
    int dynamic_instability_flag = 0;
    double diameter = 1;
    double length = -1;
    double persistence_length = 400;
    double max_length = 500;
    double min_length = 1;
    double max_bond_length = 5;
    double min_bond_length = 1.5;
    int spiral_flag = 0;
    double driving_factor = 0;
    double friction_ratio = 2;
    int force_induced_catastrophe_flag = 0;
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
    double packing_fraction = -1;
    int shuffle = 1;
    double shuffle_factor = 10;
    double shuffle_frequency = 0.001;
    double perlen_ratio = -1;
    int n_bonds = -1;
    int driving_method = 0;
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
    int n_filaments_mother = 0;
    double diameter = 10;
    double length = 20;
    int n_filaments_bud = 1;
    int alignment_potential = 0;
    double k_spring = 1000;
    double k_align = 0;
    double spring_length = 0;
    double spb_diameter = 5;
};

class motor_parameters : public species_parameters {
  public:
    double k_off = 2;
    double k_on = 10;
    double concentration = 0;
    double velocity = 1;
    int diffusion_flag = 1;
    double k_spring = 1000;
    double f_spring_max = 100;
    double diameter = 3;
    int walker = 0;
    int step_direction = 0;
};

class system_parameters {
  public:
    long seed = 7859459105545;
    int n_runs = 1;
    int n_random = 1;
    std::string run_name = "sc";
    int n_dim = 3;
    int n_periodic = 0;
    int boundary = 0;
    double system_radius = 100;
    double delta = 0.001;
    int graph_background = 1;
    long n_steps = 1000000;
    double cell_length = 10;
    int n_update_cells = 0;
    int graph_flag = 0;
    int n_graph = 1000;
    double graph_diameter = 0;
    int load_checkpoint = 0;
    int movie_flag = 0;
    int draw_boundary = 1;
    std::string checkpoint_run_name = "sc";
    std::string insertion_type = "species";
    int time_flag = 0;
    std::string movie_directory = "frames";
    double bud_height = 680;
    double bud_radius = 300;
    double lj_epsilon = 1;
    double wca_eps = 1;
    double wca_sig = 1;
    double target_pressure = 0;
    double ss_a = 1;
    double ss_rs = 1.5;
    double ss_eps = 1;
    double f_cutoff = 100;
    int max_overlap = 100000;
    int constant_volume = 0;
    int constant_pressure = 0;
    int pressure_time = 100;
    double target_radius = 100;
    double compressibility = 1;
    int stoch_flag = 1;
    int thermo_flag = 0;
    int n_thermo = 1000;
    int interaction_flag = 1;
    int species_insertion_failure_threshold = 10000;
    int uniform_crystal = 0;
    int n_steps_equil = 0;
    int static_particle_number = 0;
    species_parameters species;
    filament_parameters filament;
    hard_rod_parameters hard_rod;
    br_bead_parameters br_bead;
    md_bead_parameters md_bead;
    centrosome_parameters centrosome;
    bead_spring_parameters bead_spring;
    spherocylinder_parameters spherocylinder;
    spindle_parameters spindle;
    motor_parameters motor;
};

#endif // _SIMCORE_PARAMETERS_H_