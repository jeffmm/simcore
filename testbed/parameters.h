#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_

struct system_parameters {
  unsigned int n_dim;
  double system_radius;
  double delta;
  boundary_type boundary;
  int graph_background;
  int n_steps;
};

#endif
