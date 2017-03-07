#ifndef _SIMCORE_PARSE_FLAGS_H_
#define _SIMCORE_PARSE_FLAGS_H_
#include <getopt.h>
#include <string>
#include <vector>

// Run parameters that needs to be carried forward should be stored here.
struct run_options {
  int n_runs = 1;
  int debug = 0;
  int r_flag = 0;
  int n_flag = 0;
  int g_flag = 0;
  int m_flag = 0;
  int a_flag = 0;
  int p_flag = 0;
  int l_flag = 0;
  int n_graph = 100;
  std::string param_file;
  std::string run_name = "sc";
};

// Define flags here
// NOTE: Need to update n_flags, PARSE_OPTS function case-switches, 
//       as well as any flag descriptions in the static string array below
//       when adding new flags.

static const int n_flags = 9;
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"debug", no_argument, 0, 'd'},
  {"run-name", required_argument, 0, 'r'},
  {"n-runs", required_argument, 0, 'n'},
  {"graphics", required_argument, 0, 'g'},
  {"movie", no_argument, 0, 'm'},
  {"analysis", no_argument, 0, 'a'},
  {"posit", no_argument, 0, 'p'},
  {"load",no_argument, 0, 'l'},
  {0, 0, 0, 0}
};

// Descriptions for flags
static const std::string desc[n_flags][2] = {
  {"show this menu\n", "none"},
  {"run program in debug mode\n", "none"},
  {"where rname is the name of a run session (for organizing batch jobs)\n", "rname"},
  {"where num is the number of independent runs to perform with the given parameters\n", "num"},
  {"run graphics from posit/spec files without recording, with n_graph = ngraph\n", "ngraph"},
  {"run simulation to make movie using spec files with same run-name\n", "none"},
  {"run analysis on spec files with same run-name for parameter species", "none"},
  {"use posit files for movies/analysis rather than spec files\n", "none"},
  {"run simulation starting from checkpoint files with same run-name for parameter species\n", "none"}
};


/*************************
   SHOW_HELP_INFO
    Prints help info when an invalid command is given or when the
    help flag is given.
**************************/
static void show_help_info(std::string progname) {
  std::cout << "\n" << "  SimCORE usage: " << std::endl;
  std::cout << "    " << progname << " params_file --flag1 option1 --flag2 option2 ..." << std::endl;
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
static run_options parse_opts(int argc, char *argv[]) {
  
  if (argc == 1) {
    show_help_info(argv[0]);
    exit(1);
  }

  run_options run_opts;
  int tmp;
  while (1) {
    int option_index = 0;
    tmp = getopt_long(argc, argv, "hdmaplg:r:n:", long_options, &option_index);
    if (tmp == -1)
      break;
    switch (tmp) {
      case 0:
        break;
      case 'h':
        show_help_info(argv[0]);
        exit(0);
      case 'd':
        std::cout << "  Running in debug mode\n";
        run_opts.debug = 1;
        break;
      case 'r':
        run_opts.r_flag = 1;
        run_opts.run_name = optarg;
        break;
      case 'n':
        run_opts.n_flag = 1;
        run_opts.n_runs = atoi(optarg);
        break;
      case 'g':
        run_opts.g_flag = 1;
        run_opts.n_graph = atoi(optarg);
        break;
      case 'm':
        run_opts.m_flag = 1;
        break;
      case 'a':
        run_opts.a_flag = 1;
        break;
      case 'p':
        run_opts.p_flag = 1;
        break;
      case 'l':
        run_opts.l_flag = 1;
        break;
      case '?':
        exit(1);
      default:
        std::cout << "  ERROR: Unrecognized option: " << tmp << "\n";
        exit(1);
    }
  }
  if (optind < argc) {
    run_opts.param_file = argv[optind];
    printf("  Using parameter file: %s\n",argv[optind++]);
  }
  else {
    printf("ERROR! No parameter file given! Exiting.\n");
    exit(1);
  }
  if (optind < argc) {
    printf("ERROR: Unrecognized context for non-option arguments: ");
    while (optind<argc) 
      printf("%s ",argv[optind++]);
    putchar('\n');
    exit(1);
  }
  if (run_opts.m_flag) {
    printf("  Making movies");
    if (run_opts.p_flag) {
      printf(" using posit files.");
    }
  }
  putchar('\n');
  if (run_opts.a_flag) {
    printf("  Running analyses");
    if (run_opts.p_flag) {
      printf(" using posit files.");
    }
  }
  putchar('\n');
    //if (run_opts.a_flag!=1 && run_opts.m_flag!=1) {
      //printf("ERROR: Unrecognized context for non-option arguments: ");
      //while (optind < argc)
        //printf ("%s ", argv[optind++]);
      //putchar ('\n');
      //exit(1);
    //}
    //else {
      //printf("Posit files: ");
      //while (optind < argc) {
        //run_opts.posit_files.push_back(argv[optind++]);
      //}
      //for (int i=0; i<run_opts.posit_files.size(); ++i)
        //std::cout << run_opts.posit_files[i] << " ";
      //putchar('\n');
    //}
  //}
  return run_opts;
}

#endif // _SIMCORE_PARSE_FLAGS_H_
