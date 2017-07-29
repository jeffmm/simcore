#ifndef _SIMCORE_PARAMETERS_H_
#define _SIMCORE_PARAMETERS_H_
#include "definitions.h"

struct system_parameters {
  unsigned int n_dim;
  double system_radius;
  boundary_type boundary;
  int graph_background;
};

#endif
