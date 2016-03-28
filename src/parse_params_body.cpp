if ( param_name.compare("seed") == 0 ) {
  params_.seed = atol(param_value.c_str()); 
  std::cout << "  seed = " << param_value << std::endl;
}
if ( param_name.compare("n_dim") == 0 ) {
  params_.n_dim = atoi(param_value.c_str()); 
  std::cout << "  n_dim = " << param_value << std::endl;
}
if ( param_name.compare("n_periodic") == 0 ) {
  params_.n_periodic = atoi(param_value.c_str()); 
  std::cout << "  n_periodic = " << param_value << std::endl;
}
if ( param_name.compare("insert_type") == 0 ) {
  params_.insert_type = atoi(param_value.c_str()); 
  std::cout << "  insert_type = " << param_value << std::endl;
}
if ( param_name.compare("boundary_type") == 0 ) {
  params_.boundary_type = atoi(param_value.c_str()); 
  std::cout << "  boundary_type = " << param_value << std::endl;
}
if ( param_name.compare("n_filaments_free") == 0 ) {
  params_.n_filaments_free = atoi(param_value.c_str()); 
  std::cout << "  n_filaments_free = " << param_value << std::endl;
}
if ( param_name.compare("n_filaments_attached") == 0 ) {
  params_.n_filaments_attached = atoi(param_value.c_str()); 
  std::cout << "  n_filaments_attached = " << param_value << std::endl;
}
if ( param_name.compare("n_spheres") == 0 ) {
  params_.n_spheres = atoi(param_value.c_str()); 
  std::cout << "  n_spheres = " << param_value << std::endl;
}
if ( param_name.compare("n_steps") == 0 ) {
  params_.n_steps = atoi(param_value.c_str()); 
  std::cout << "  n_steps = " << param_value << std::endl;
}
if ( param_name.compare("n_graph") == 0 ) {
  params_.n_graph = atoi(param_value.c_str()); 
  std::cout << "  n_graph = " << param_value << std::endl;
}
if ( param_name.compare("dynamic_instability_flag") == 0 ) {
  params_.dynamic_instability_flag = atoi(param_value.c_str()); 
  std::cout << "  dynamic_instability_flag = " << param_value << std::endl;
}
if ( param_name.compare("force_induced_catastrophe_flag") == 0 ) {
  params_.force_induced_catastrophe_flag = atoi(param_value.c_str()); 
  std::cout << "  force_induced_catastrophe_flag = " << param_value << std::endl;
}
if ( param_name.compare("error_analysis_flag") == 0 ) {
  params_.error_analysis_flag = atoi(param_value.c_str()); 
  std::cout << "  error_analysis_flag = " << param_value << std::endl;
}
if ( param_name.compare("graph_flag") == 0 ) {
  params_.graph_flag = atoi(param_value.c_str()); 
  std::cout << "  graph_flag = " << param_value << std::endl;
}
if ( param_name.compare("grab_flag") == 0 ) {
  params_.grab_flag = atoi(param_value.c_str()); 
  std::cout << "  grab_flag = " << param_value << std::endl;
}
if ( param_name.compare("pair_interaction_flag") == 0 ) {
  params_.pair_interaction_flag = atoi(param_value.c_str()); 
  std::cout << "  pair_interaction_flag = " << param_value << std::endl;
}
if ( param_name.compare("buckling_analysis_flag") == 0 ) {
  params_.buckling_analysis_flag = atoi(param_value.c_str()); 
  std::cout << "  buckling_analysis_flag = " << param_value << std::endl;
}
if ( param_name.compare("time_analysis_flag") == 0 ) {
  params_.time_analysis_flag = atoi(param_value.c_str()); 
  std::cout << "  time_analysis_flag = " << param_value << std::endl;
}
if ( param_name.compare("theta_validation_flag") == 0 ) {
  params_.theta_validation_flag = atoi(param_value.c_str()); 
  std::cout << "  theta_validation_flag = " << param_value << std::endl;
}
if ( param_name.compare("save_state_flag") == 0 ) {
  params_.save_state_flag = atoi(param_value.c_str()); 
  std::cout << "  save_state_flag = " << param_value << std::endl;
}
if ( param_name.compare("n_save_state") == 0 ) {
  params_.n_save_state = atoi(param_value.c_str()); 
  std::cout << "  n_save_state = " << param_value << std::endl;
}
if ( param_name.compare("n_bins") == 0 ) {
  params_.n_bins = atoi(param_value.c_str()); 
  std::cout << "  n_bins = " << param_value << std::endl;
}
if ( param_name.compare("n_validate") == 0 ) {
  params_.n_validate = atoi(param_value.c_str()); 
  std::cout << "  n_validate = " << param_value << std::endl;
}
if ( param_name.compare("n_buckle_on") == 0 ) {
  params_.n_buckle_on = atoi(param_value.c_str()); 
  std::cout << "  n_buckle_on = " << param_value << std::endl;
}
if ( param_name.compare("n_buckle_off") == 0 ) {
  params_.n_buckle_off = atoi(param_value.c_str()); 
  std::cout << "  n_buckle_off = " << param_value << std::endl;
}
if ( param_name.compare("metric_forces") == 0 ) {
  params_.metric_forces = atoi(param_value.c_str()); 
  std::cout << "  metric_forces = " << param_value << std::endl;
}
if ( param_name.compare("graph_background") == 0 ) {
  params_.graph_background = atoi(param_value.c_str()); 
  std::cout << "  graph_background = " << param_value << std::endl;
}

if ( param_name.compare("delta") == 0 ) {
  params_.delta = atof(param_value.c_str()); 
  std::cout << "  delta = " << param_value << std::endl;
}
if ( param_name.compare("friction_ratio") == 0 ) {
  params_.friction_ratio = atof(param_value.c_str()); 
  std::cout << "  friction_ratio = " << param_value << std::endl;
}
if ( param_name.compare("mother_daughter_dist") == 0 ) {
  params_.mother_daughter_dist = atof(param_value.c_str()); 
  std::cout << "  mother_daughter_dist = " << param_value << std::endl;
}
if ( param_name.compare("mother_radius") == 0 ) {
  params_.mother_radius = atof(param_value.c_str()); 
  std::cout << "  mother_radius = " << param_value << std::endl;
}
if ( param_name.compare("daughter_radius") == 0 ) {
  params_.daughter_radius = atof(param_value.c_str()); 
  std::cout << "  daughter_radius = " << param_value << std::endl;
}
if ( param_name.compare("sphere_radius") == 0 ) {
  params_.sphere_radius = atof(param_value.c_str()); 
  std::cout << "  sphere_radius = " << param_value << std::endl;
}
if ( param_name.compare("filament_diameter") == 0 ) {
  params_.filament_diameter = atof(param_value.c_str()); 
  std::cout << "  filament_diameter = " << param_value << std::endl;
}
if ( param_name.compare("max_length") == 0 ) {
  params_.max_length = atof(param_value.c_str()); 
  std::cout << "  max_length = " << param_value << std::endl;
}
if ( param_name.compare("min_length") == 0 ) {
  params_.min_length = atof(param_value.c_str()); 
  std::cout << "  min_length = " << param_value << std::endl;
}
if ( param_name.compare("max_segment_length") == 0 ) {
  params_.max_segment_length = atof(param_value.c_str()); 
  std::cout << "  max_segment_length = " << param_value << std::endl;
}
if ( param_name.compare("min_segment_length") == 0 ) {
  params_.min_segment_length = atof(param_value.c_str()); 
  std::cout << "  min_segment_length = " << param_value << std::endl;
}
if ( param_name.compare("persistence_length") == 0 ) {
  params_.persistence_length = atof(param_value.c_str()); 
  std::cout << "  persistence_length = " << param_value << std::endl;
}
if ( param_name.compare("spring_filament_sphere") == 0 ) {
  params_.spring_filament_sphere = atof(param_value.c_str()); 
  std::cout << "  spring_filament_sphere = " << param_value << std::endl;
}
if ( param_name.compare("spring_buckling_init") == 0 ) {
  params_.spring_buckling_init = atof(param_value.c_str()); 
  std::cout << "  spring_buckling_init = " << param_value << std::endl;
}
if ( param_name.compare("buckle_rate") == 0 ) {
  params_.buckle_rate = atof(param_value.c_str()); 
  std::cout << "  buckle_rate = " << param_value << std::endl;
}
if ( param_name.compare("buckle_parameter") == 0 ) {
  params_.buckle_parameter = atof(param_value.c_str()); 
  std::cout << "  buckle_parameter = " << param_value << std::endl;
}
if ( param_name.compare("r_cutoff_sphere") == 0 ) {
  params_.r_cutoff_sphere = atof(param_value.c_str()); 
  std::cout << "  r_cutoff_sphere = " << param_value << std::endl;
}
if ( param_name.compare("r_cutoff_boundary") == 0 ) {
  params_.r_cutoff_boundary = atof(param_value.c_str()); 
  std::cout << "  r_cutoff_boundary = " << param_value << std::endl;
}
if ( param_name.compare("f_shrink_to_grow") == 0 ) {
  params_.f_shrink_to_grow = atof(param_value.c_str()); 
  std::cout << "  f_shrink_to_grow = " << param_value << std::endl;
}
if ( param_name.compare("f_pause_to_grow") == 0 ) {
  params_.f_pause_to_grow = atof(param_value.c_str()); 
  std::cout << "  f_pause_to_grow = " << param_value << std::endl;
}
if ( param_name.compare("f_grow_to_shrink") == 0 ) {
  params_.f_grow_to_shrink = atof(param_value.c_str()); 
  std::cout << "  f_grow_to_shrink = " << param_value << std::endl;
}
if ( param_name.compare("f_grow_to_pause") == 0 ) {
  params_.f_grow_to_pause = atof(param_value.c_str()); 
  std::cout << "  f_grow_to_pause = " << param_value << std::endl;
}
if ( param_name.compare("f_pause_to_shrink") == 0 ) {
  params_.f_pause_to_shrink = atof(param_value.c_str()); 
  std::cout << "  f_pause_to_shrink = " << param_value << std::endl;
}
if ( param_name.compare("f_shrink_to_pause") == 0 ) {
  params_.f_shrink_to_pause = atof(param_value.c_str()); 
  std::cout << "  f_shrink_to_pause = " << param_value << std::endl;
}
if ( param_name.compare("v_poly") == 0 ) {
  params_.v_poly = atof(param_value.c_str()); 
  std::cout << "  v_poly = " << param_value << std::endl;
}
if ( param_name.compare("v_depoly") == 0 ) {
  params_.v_depoly = atof(param_value.c_str()); 
  std::cout << "  v_depoly = " << param_value << std::endl;
}
if ( param_name.compare("graph_diameter") == 0 ) {
  params_.graph_diameter = atof(param_value.c_str()); 
  std::cout << "  graph_diameter = " << param_value << std::endl;
}

if ( param_name.compare("grab_file") == 0 ) {
  params_.grab_file = (char *) gmalloc((strlen(param_value.c_str()) + 1) * sizeof(char));
  strcpy(params_.grab_file, param_value.c_str());
  std::cout << "  grab_file = " << param_value << std::endl;
}

if ( param_name.compare("cell_list_flag") == 0 ) {
  params_.cell_list_flag = atoi(param_value.c_str()); 
  std::cout << "  cell_list_flag = " << param_value << std::endl;
}

if ( param_name.compare("potential_flag") == 0){
    params_.potential_flag = atoi(param_value.c_str());
    std::cout << "  potential_flag = " << param_value << std::endl;
}

if ( param_name.compare("validate_file_number") == 0){
    params_.validate_file_number = atoi(param_value.c_str());
    std::cout << "  validate_file_number = " << param_value << std::endl;
}

if (param_name.compare("n_particles") == 0){
    params_.n_particles = atoi(param_value.c_str());
    std::cout << "  n_particles = " << param_value << std::endl;
}

if (param_name.compare("particle_radius") == 0){
    params_.particle_radius = atof(param_value.c_str());
    std::cout << "  particle_radius = " << param_value << std::endl;
}

if (param_name.compare("temp") == 0){
    params_.temp = atof(param_value.c_str());
    std::cout << "  temperature = " << param_value << std::endl;
}
