#ifndef _CYTOSCORE_PARAMETERS_H_
#define _CYTOSCORE_PARAMETERS_H_

struct system_parameters { 

  long seed; // seed for random number generator
  int n_dim,
      n_periodic,
      insert_type, // 0 for Random, manual otherwise
      boundary_type,
      n_filaments_free, // number of free filaments in system
      n_filaments_attached, // number of attached filaments (per sphere)
      n_spheres, // number of spheres in system
      n_steps, // number of steps for simulation
      n_graph, // number of frames between refreshing graphics
      dynamic_instability_flag,
      force_induced_catastrophe_flag,
      error_analysis_flag,
      graph_flag,
      grab_flag,
      pair_interaction_flag,
      buckling_analysis_flag,
      time_analysis_flag,
      theta_validation_flag,
      save_state_flag,
      cell_list_flag,
      potential_flag,
      n_save_state, // number of frames between saving the current system state, useful for debugging
      n_bins,
      n_validate,
      n_buckle_on,
      n_buckle_off,
      metric_forces,
      validate_file_number,
      graph_background,
      n_particles;
  double delta, // time duration between each simulation step
         friction_ratio,
         mother_daughter_dist, // distance between mother-daughter cell centers
         mother_radius, // radius of mother cell
         daughter_radius, // radius of daughter cell
         sphere_radius, // radius of sphere
         filament_diameter, // diameter of filament
         min_length,
         max_length,
         min_segment_length,
         max_segment_length,
         persistence_length, // persistence_length of microtubules
         spring_filament_sphere, // spring constant for filament sphere attachments
         spring_buckling_init, // initial spring constant for buckling forces
         buckle_rate,
         buckle_parameter,
         r_cutoff_sphere, // cutoff for sphere wca
         r_cutoff_boundary, // cutoff for boundary wca
         f_shrink_to_grow, // Frequencies for dynamic instability
         f_shrink_to_pause,
         f_pause_to_grow,
         f_pause_to_shrink,
         f_grow_to_shrink,
         f_grow_to_pause,
         v_poly,
         v_depoly,
         graph_diameter,
         particle_radius,
         particle_mass,
         temp;
  
  char *grab_file; // root directory to save bmp image of graphics to

// Initialization of parameters to default values 
  inline void init() {
#include "parameter_defaults.h"
  }

};


#endif // _CYTOSCORE_PARAMETERS_H_
