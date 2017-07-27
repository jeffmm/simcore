#include "object.h"

int main() {

  unsigned int n_dim = 3;
  Object::SetNDim(n_dim);
  Mesh m;
  double pos[3] = {0,0,0};
  double u[3] = {0,0,1};
  m.InitBondAt(pos,u,2,1);
  m.Report();
  return 0;
}

