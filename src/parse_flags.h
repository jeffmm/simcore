#ifndef _SIMCORE_PARSE_FLAGS_H_
#define _SIMCORE_PARSE_FLAGS_H_
#include <getopt.h>
#include <string>

// Define flags here
// NOTE: Need to update n_flags, PARSE_OPTS function case-switches, 
//       as well as any flag descriptions in the static string array below
//       when adding new flags.

static const int n_flags = 7; 
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"debug", no_argument, 0, 'd'},
  {"file", required_argument, 0, 'f'},
  {"run-name", required_argument, 0, 'r'},
  {"n-runs", required_argument, 0, 'n'},
  {"test", no_argument, 0, 't'},
  {"movie", required_argument, 0, 'm'},
  {0, 0, 0, 0}
};

// Descriptions for flags
static const std::string desc[n_flags][2] = {
  {"show this menu\n", "none"},
  {"run program in debug mode\n", "none"},
  {"where fname is the input parameter file (REQUIRED)\n", "fname"},
  {"where rname is the name of a run session (for organizing batch jobs)\n", "rname"},
  {"where num is the number of independent runs to perform with the given parameters\n", "num"},
  {"where movie is the posit file name to recreate previous simulation\n", "movie"},
  {"run program in test mode\n", "none"}
};

// Run parameters that needs to be carried forward should be stored here.
struct run_options {
  run_options() {n_runs = 1; debug = 0; test = 0; f_flag = 0; r_flag = 0; n_flag = 0; m_flag = 0; run_name = "sc";}
  int n_runs;
  int debug;
  int test;
  int f_flag;
  int r_flag;
  int n_flag;
  int m_flag;
  std::string param_file;
  std::string run_name;
  std::string posit_file;
};

/*************************
   SHOW_HELP_INFO
    Prints help info when an invalid command is given or when the
    help flag is given.
**************************/
void show_help_info(std::string progname) {
  std::cout << "\n" << "  CytoSCORe usage: " << std::endl;
  std::cout << "    " << progname << " --flag1 option1 --flag2 option2 ..." << std::endl;
  std::cout << "\n" << "  where flags are one of the following: " << std::endl;
  for (int i=0; i<n_flags; ++i) {
    std::cout << "    --" << long_options[i].name;
    if (desc[i][1].compare("none") != 0)
      std::cout << " " << desc[i][1];
    std::cout << ", -" << (char) long_options[i].val;
    if (desc[i][1].compare("none") != 0)
      std::cout << " " << desc[i][1];
    std::cout << ", " << desc[i][0];
  }
}

/*************************
   PARSE_OPTS
    New flag switches should be added below,
    as well as any information that needs to be carried forward
    in the run_options structure defined above.
**************************/
run_options parse_opts(int argc, char *argv[]) {
  
  if (argc == 1) {
    show_help_info(argv[0]);
    exit(1);
  }

  run_options run_opts;
  int tmp;
  while (1) {
    int option_index = 0;
    tmp = getopt_long(argc, argv, "hdtr:f:n:m:", long_options, &option_index);
    if (tmp == -1)
      break;
    
    switch (tmp) {
      case 0:
        break;
      case 'h':
        show_help_info(argv[0]);
        exit(0);
      case 'd':
        std::cout << "  Running with Debug Mode\n";
        run_opts.debug = 1;
        break;
      case 't':
        std::cout << "  Running Test Mode\n";
        run_opts.test = 1;
        break;
      case 'f':
        run_opts.f_flag = 1;
        run_opts.param_file = optarg;
        break;
      case 'r':
        run_opts.r_flag = 1;
        run_opts.run_name = optarg;
        break;
      case 'n':
        run_opts.n_flag = 1;
        run_opts.n_runs = atoi(optarg);
        break;
      case 'm':
        run_opts.m_flag = 1;
        run_opts.posit_file = optarg;
        break;
      case '?':
        break;
      default:
        std::cout << "  ERROR: Unrecognized option!\n";
        exit(1);
    }
  }
  return run_opts;
}

#endif // _SIMCORE_PARSE_FLAGS_H_
