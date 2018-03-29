//#define NOBJS 50000 // Needs to be <= 500000 to fit on the heap
//#define NDIM 2
//#define NPER 2
//#define SEED 819844816414
//#define BOX_SIZE 900
//#define RCUT 1
//#define DEBUG 0
//#define NINT(x)            ((x) < 0.0 ? (int) ((x) - 0.5) : (int) ((x) + 0.5))
//#define FFF_NUM 100
//#define FFF_AR 100
//#define FFF_LP_RATIO 10
//#define FFF_NBONDS 25
//#define FFF_DRIVING 0
//#define FFF_STOCH 1
//#define FFF_METRIC 1
//#define FFF_FRICTION_RATIO 2
//#define FFF_INSERT 0 // 0 for isotropic, 1 for nematic

struct system_params {
  int n_dim = 2;
  int n_per = 2;
  long seed = 777777;
  double system_diameter = 900;
  double cell_length = 1;
  int debug = 1;
  int n_filaments = 100;
  double aspect_ratio = 100;
  double lp_ratio = 100;
  int n_bonds = 25;
  double driving = 0;
  int stoch = 1;
  int metric = 1;
  double friction_ratio = 2;
  int insert_type = 0;
  int n_sites;
};
