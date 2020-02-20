
double cpu_time();
void grabber(int width, int height, std::string fname, int framenum);
void rotate_vector_relative(int n_dim, double *vect1, double *vect2);
void rotate_vector(double *v, double *k, double theta, int n_dim);
double dot_product(int n_dim, double const *const a, double const *const b);
double determinant(int n, double **mat);
void separation_vector(int n_dim, int n_periodic, double const *const r1,
                       double const *const s1, double const *const r2,
                       double const *const s2, double **unit_cell, double *dr);
void cross_product(double const *const a, double const *const b, double *c,
                   int n_dim);
void normalize_vector(double *a, int n_dim);
void tridiagonal_solver(std::vector<double> *a, std::vector<double> *b,
                        std::vector<double> *c, std::vector<double> *d, int n);
void invert_sym_2d_matrix(double *a, double *b);
void invert_sym_3d_matrix(double *a, double *b);
void periodic_boundary_conditions(int n_dim, int n_periodic, double *h,
                                  double *h_inv, double *r, double *s);
