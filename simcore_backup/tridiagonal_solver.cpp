#include "auxiliary.h"

void tridiagonal_solver(std::vector<double> *a, std::vector<double> *b, std::vector<double> *c, std::vector<double> *d, int n) {
  
  //  Solves Ax=d for x as seen below:
  //
  //  | b[0] c[0]   0  | | x[0] |   | d[0] |
  //  | a[0] b[1] c[1] | | x[1] | = | d[1] |
  //  |   0  a[1] b[2] | | x[2] |   | d[2] |

  //  Works only for tridiagonal matrices.
  //  Uses Thomas' algorithm which is unstable
  //  unless matrix is guaranteed to be symmetric
  //  positive definite.

  //  Input: pointer to array of lower off-diagonal elements (*a)
  //         pointer to array of diagonal elements as seen above (*b)
  //         pointer to array of upper off-diagonal elements  (*c)
  //         pointer to array of d vector (*d) 
  //         number of unknowns (n)
  //
  //  Output: Array of solutions x (*d)

  n--;
  (*c)[0]/=(*b)[0];
  (*d)[0]/=(*b)[0];
  for (int i=1; i<n; ++i) {
    (*c)[i] /= (*b)[i] - (*a)[i-1]*(*c)[i-1];
    (*d)[i] = ((*d)[i] - (*a)[i-1]*(*d)[i-1]) / ((*b)[i] - (*a)[i-1]*(*c)[i-1]);
  }
  (*d)[n] = ((*d)[n] - (*a)[n-1]*(*d)[n-1]) / ((*b)[n] - (*a)[n-1]*(*c)[n-1]);
  for (int i=n; i-- >0;) 
    (*d)[i] -= (*c)[i]*(*d)[i+1];
  return;
}
