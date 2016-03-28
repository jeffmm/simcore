
#include <iostream>
#include <stdlib.h>
#include "parse_flags.h"
#include "simulation_manager.h"
/*************************
   MAIN
**************************/
int main(int argc, char *argv[]) {

  // Parse input flags, see parse_flags.h for documentation
  run_options run_opts = parse_opts(argc, argv);
  if (!run_opts.file_flag && !run_opts.debug) {
    std::cout << "  ERROR: No parameter file given!\n";
    show_help_info(argv[0]);
    exit(1);
  }
  if (run_opts.debug) {
    SimulationManager sim;
  }
  else {
    SimulationManager sim(run_opts.param_file, run_opts.run_name, run_opts.n_runs);
  }
  return 0;
}

