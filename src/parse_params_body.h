// parse_params_body.h, generated automatically using make_params

if (param_name.compare("n_runs") == 0 || param_name.compare("run_name") == 0) {}
else if (param_name.compare("seed") == 0){
  params_[i_var].seed = atol(param_value.c_str());
}
else if (param_name.compare("n_dim") == 0){
  params_[i_var].n_dim = atoi(param_value.c_str());
}
else if (param_name.compare("n_periodic") == 0){
  params_[i_var].n_periodic = atoi(param_value.c_str());
}
else if (param_name.compare("boundary_type") == 0){
  params_[i_var].boundary_type = atoi(param_value.c_str());
}
else if (param_name.compare("system_radius") == 0){
  params_[i_var].system_radius = atof(param_value.c_str());
}
else if (param_name.compare("insert_type") == 0){
  params_[i_var].insert_type = atoi(param_value.c_str());
}
else if (param_name.compare("n_filaments_free") == 0){
  params_[i_var].n_filaments_free = atoi(param_value.c_str());
}
else if (param_name.compare("n_filaments_attached") == 0){
  params_[i_var].n_filaments_attached = atoi(param_value.c_str());
}
else if (param_name.compare("n_spheres") == 0){
  params_[i_var].n_spheres = atoi(param_value.c_str());
}
else if (param_name.compare("n_particles") == 0){
  params_[i_var].n_particles = atoi(param_value.c_str());
}
else if (param_name.compare("particle_radius") == 0){
  params_[i_var].particle_radius = atof(param_value.c_str());
}
else if (param_name.compare("particle_mass") == 0){
  params_[i_var].particle_mass = atof(param_value.c_str());
}
else if (param_name.compare("n_steps") == 0){
  params_[i_var].n_steps = atoi(param_value.c_str());
}
else if (param_name.compare("n_graph") == 0){
  params_[i_var].n_graph = atoi(param_value.c_str());
}
else if (param_name.compare("dynamic_instability_flag") == 0){
  params_[i_var].dynamic_instability_flag = atoi(param_value.c_str());
}
else if (param_name.compare("mother_daughter_dist") == 0){
  params_[i_var].mother_daughter_dist = atof(param_value.c_str());
}
else if (param_name.compare("force_induced_catastrophe_flag") == 0){
  params_[i_var].force_induced_catastrophe_flag = atoi(param_value.c_str());
}
else if (param_name.compare("error_analysis_flag") == 0){
  params_[i_var].error_analysis_flag = atoi(param_value.c_str());
}
else if (param_name.compare("graph_flag") == 0){
  params_[i_var].graph_flag = atoi(param_value.c_str());
}
else if (param_name.compare("grab_flag") == 0){
  params_[i_var].grab_flag = atoi(param_value.c_str());
}
else if (param_name.compare("pair_interaction_flag") == 0){
  params_[i_var].pair_interaction_flag = atoi(param_value.c_str());
}
else if (param_name.compare("theta_validation_flag") == 0){
  params_[i_var].theta_validation_flag = atoi(param_value.c_str());
}
else if (param_name.compare("buckling_analysis_flag") == 0){
  params_[i_var].buckling_analysis_flag = atoi(param_value.c_str());
}
else if (param_name.compare("time_flag") == 0){
  params_[i_var].time_flag = atoi(param_value.c_str());
}
else if (param_name.compare("save_state_flag") == 0){
  params_[i_var].save_state_flag = atoi(param_value.c_str());
}
else if (param_name.compare("n_save_state") == 0){
  params_[i_var].n_save_state = atoi(param_value.c_str());
}
else if (param_name.compare("n_bins") == 0){
  params_[i_var].n_bins = atoi(param_value.c_str());
}
else if (param_name.compare("n_validate") == 0){
  params_[i_var].n_validate = atoi(param_value.c_str());
}
else if (param_name.compare("n_buckle_on") == 0){
  params_[i_var].n_buckle_on = atoi(param_value.c_str());
}
else if (param_name.compare("n_buckle_off") == 0){
  params_[i_var].n_buckle_off = atoi(param_value.c_str());
}
else if (param_name.compare("metric_forces") == 0){
  params_[i_var].metric_forces = atoi(param_value.c_str());
}
else if (param_name.compare("delta") == 0){
  params_[i_var].delta = atof(param_value.c_str());
}
else if (param_name.compare("friction_ratio") == 0){
  params_[i_var].friction_ratio = atof(param_value.c_str());
}
else if (param_name.compare("daughter_radius") == 0){
  params_[i_var].daughter_radius = atof(param_value.c_str());
}
else if (param_name.compare("sphere_radius") == 0){
  params_[i_var].sphere_radius = atof(param_value.c_str());
}
else if (param_name.compare("filament_diameter") == 0){
  params_[i_var].filament_diameter = atof(param_value.c_str());
}
else if (param_name.compare("min_length") == 0){
  params_[i_var].min_length = atof(param_value.c_str());
}
else if (param_name.compare("max_length") == 0){
  params_[i_var].max_length = atof(param_value.c_str());
}
else if (param_name.compare("min_segment_length") == 0){
  params_[i_var].min_segment_length = atof(param_value.c_str());
}
else if (param_name.compare("max_segment_length") == 0){
  params_[i_var].max_segment_length = atof(param_value.c_str());
}
else if (param_name.compare("persistence_length") == 0){
  params_[i_var].persistence_length = atof(param_value.c_str());
}
else if (param_name.compare("spring_filament_sphere") == 0){
  params_[i_var].spring_filament_sphere = atof(param_value.c_str());
}
else if (param_name.compare("spring_buckling_init") == 0){
  params_[i_var].spring_buckling_init = atof(param_value.c_str());
}
else if (param_name.compare("buckle_rate") == 0){
  params_[i_var].buckle_rate = atof(param_value.c_str());
}
else if (param_name.compare("buckle_parameter") == 0){
  params_[i_var].buckle_parameter = atof(param_value.c_str());
}
else if (param_name.compare("r_cutoff_sphere") == 0){
  params_[i_var].r_cutoff_sphere = atof(param_value.c_str());
}
else if (param_name.compare("r_cutoff_boundary") == 0){
  params_[i_var].r_cutoff_boundary = atof(param_value.c_str());
}
else if (param_name.compare("f_shrink_to_grow") == 0){
  params_[i_var].f_shrink_to_grow = atof(param_value.c_str());
}
else if (param_name.compare("f_shrink_to_pause") == 0){
  params_[i_var].f_shrink_to_pause = atof(param_value.c_str());
}
else if (param_name.compare("f_pause_to_grow") == 0){
  params_[i_var].f_pause_to_grow = atof(param_value.c_str());
}
else if (param_name.compare("f_pause_to_shrink") == 0){
  params_[i_var].f_pause_to_shrink = atof(param_value.c_str());
}
else if (param_name.compare("f_grow_to_pause") == 0){
  params_[i_var].f_grow_to_pause = atof(param_value.c_str());
}
else if (param_name.compare("f_grow_to_shrink") == 0){
  params_[i_var].f_grow_to_shrink = atof(param_value.c_str());
}
else if (param_name.compare("v_poly") == 0){
  params_[i_var].v_poly = atof(param_value.c_str());
}
else if (param_name.compare("v_depoly") == 0){
  params_[i_var].v_depoly = atof(param_value.c_str());
}
else if (param_name.compare("graph_diameter") == 0){
  params_[i_var].graph_diameter = atof(param_value.c_str());
}
else if ( param_name.compare("grab_file") == 0 ) {
  params_[i_var].grab_file = (char *) gmalloc((strlen(param_value.c_str()) + 1) * sizeof(char));
  strcpy(params_[i_var].grab_file, param_value.c_str());
}
else if (param_name.compare("cell_list_flag") == 0){
  params_[i_var].cell_list_flag = atoi(param_value.c_str());
}
else if (param_name.compare("rigid_tether_flag") == 0){
  params_[i_var].rigid_tether_flag = atoi(param_value.c_str());
}
else if (param_name.compare("position_correlation_flag") == 0){
  params_[i_var].position_correlation_flag = atoi(param_value.c_str());
}
else if (param_name.compare("graph_background") == 0){
  params_[i_var].graph_background = atoi(param_value.c_str());
}
else if (param_name.compare("dimer_k_spring") == 0){
  params_[i_var].dimer_k_spring = atof(param_value.c_str());
}
else if (param_name.compare("dimer_eq_length") == 0){
  params_[i_var].dimer_eq_length = atof(param_value.c_str());
}
else if (param_name.compare("dimer_length") == 0){
  params_[i_var].dimer_length = atof(param_value.c_str());
}
else if (param_name.compare("dimer_diameter") == 0){
  params_[i_var].dimer_diameter = atof(param_value.c_str());
}
else if (param_name.compare("n_dimer") == 0){
  params_[i_var].n_dimer = atoi(param_value.c_str());
}
else if (param_name.compare("n_br_bead") == 0){
  params_[i_var].n_br_bead = atoi(param_value.c_str());
}
else if (param_name.compare("br_bead_diameter") == 0){
  params_[i_var].br_bead_diameter = atof(param_value.c_str());
}
else if (param_name.compare("n_md_bead") == 0){
  params_[i_var].n_md_bead = atoi(param_value.c_str());
}
else if (param_name.compare("md_bead_diameter") == 0){
  params_[i_var].md_bead_diameter = atof(param_value.c_str());
}
else if (param_name.compare("md_bead_mass") == 0){
  params_[i_var].md_bead_mass = atof(param_value.c_str());
}
else if (param_name.compare("n_md_kmc_bead") == 0){
  params_[i_var].n_md_kmc_bead = atoi(param_value.c_str());
}
else if (param_name.compare("md_kmc_bead_diameter") == 0){
  params_[i_var].md_kmc_bead_diameter = atof(param_value.c_str());
}
else if (param_name.compare("md_kmc_bead_mass") == 0){
  params_[i_var].md_kmc_bead_mass = atof(param_value.c_str());
}
else if (param_name.compare("cell_length") == 0){
  params_[i_var].cell_length = atof(param_value.c_str());
}
else if (param_name.compare("n_update_cells") == 0){
  params_[i_var].n_update_cells = atoi(param_value.c_str());
}
else if (param_name.compare("energy_analysis_flag") == 0){
  params_[i_var].energy_analysis_flag = atoi(param_value.c_str());
}
else if (param_name.compare("lj_epsilon") == 0){
  params_[i_var].lj_epsilon = atof(param_value.c_str());
}
else if (param_name.compare("argon_diameter") == 0){
  params_[i_var].argon_diameter = atof(param_value.c_str());
}
else if (param_name.compare("neon_diameter") == 0){
  params_[i_var].neon_diameter = atof(param_value.c_str());
}
else if (param_name.compare("argon_mass") == 0){
  params_[i_var].argon_mass = atof(param_value.c_str());
}
else if (param_name.compare("neon_mass") == 0){
  params_[i_var].neon_mass = atof(param_value.c_str());
}
else if (param_name.compare("n_argon") == 0){
  params_[i_var].n_argon = atoi(param_value.c_str());
}
else if (param_name.compare("n_neon") == 0){
  params_[i_var].n_neon = atoi(param_value.c_str());
}
else if (param_name.compare("argon_rcutoff") == 0){
  params_[i_var].argon_rcutoff = atof(param_value.c_str());
}
else if (param_name.compare("neon_rcutoff") == 0){
  params_[i_var].neon_rcutoff = atof(param_value.c_str());
}
else if (param_name.compare("n_rod") == 0){
  params_[i_var].n_rod = atoi(param_value.c_str());
}
else if (param_name.compare("rod_length") == 0){
  params_[i_var].rod_length = atof(param_value.c_str());
}
else if (param_name.compare("rod_diameter") == 0){
  params_[i_var].rod_diameter = atof(param_value.c_str());
}
else if (param_name.compare("draw_interactions") == 0){
  params_[i_var].draw_interactions = atoi(param_value.c_str());
}
else if (param_name.compare("ftype") == 0){
  params_[i_var].ftype = atoi(param_value.c_str());
}
else if (param_name.compare("masterskin") == 0){
  params_[i_var].masterskin = atof(param_value.c_str());
}
else if ( param_name.compare("potfile") == 0 ) {
  params_[i_var].potfile = (char *) gmalloc((strlen(param_value.c_str()) + 1) * sizeof(char));
  strcpy(params_[i_var].potfile, param_value.c_str());
}
else {
  std::cout << "  WARNING: Parameter " << param_name << " not recognized!\n";
}
