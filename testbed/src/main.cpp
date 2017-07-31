#include "object.h"
#include "graphics.h"

int main() {

  Graphics graphics;
  //if (argc != 3) {
    //std::cerr << "ERROR! Wrong number of arguments.\n";
    //exit(1);
  //}
  int n_bonds = 100;
  long seed = 43916810496143;
  int n_dim = 3;
  Object::SetNDim(n_dim);
  Object::SetSeed(seed);
  system_parameters params;
  params.n_dim = n_dim;
  params.system_radius = 30;
  params.boundary = BOX;
  params.graph_background = 1;
  Object::SetParameters(&params);
  Mesh m;
  double pos[3] = {0,0,0};
  double u[3] = {1,0,0};
  //m.InitBondAt(pos,u,bl,1);
  m.InitSiteAt(pos,1);
  RNG rng;
  rng.Init(seed);
  for (int i=0; i<n_bonds; ++i) {
    double bl = 2+20*gsl_rng_uniform_pos(rng.r);
    m.AddRandomBondAnywhere(bl);
  }
  double background_color = (params.graph_background == 0 ? 0.1 : 1);
  std::vector<graph_struct*> graph_array;
  m.Draw(&graph_array);
  // Init space
  space_struct space;
  space.n_dim = params.n_dim;
  space.n_periodic = 0;
  space.radius = 3.0*params.system_radius;
  space.type = "box";
  double unit_cell[9] = {space.radius,0,0,0,space.radius,0,0,0,space.radius};
  //double unit_cell[4] = {space.radius,0,0,space.radius};
  space.unit_cell = unit_cell;

  graphics.Init(&graph_array, &space, background_color, 0);
  graphics.DrawLoop();

  int time=0;
  while (time<10) {
    time++;

    // ...
    // Update simulation
    // ...

    // update graphics
    graph_array.clear();
    m.Draw(&graph_array);
    graphics.Draw();
  }

  // clean up (for when we actually terminate simulation)
  graphics.Clear();

  return 0;
}

