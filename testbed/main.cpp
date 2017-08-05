#include "mesh.h"
#include "motor.h"
#include "graphics.h"

int main() {

  Graphics graphics;
  //if (argc != 3) {
    //std::cerr << "ERROR! Wrong number of arguments.\n";
    //exit(1);
  //}
  int n_bonds = 40;
  int n_motors = 20;
  long seed = 43916810496143;
  int n_dim = 3;
  Object::SetNDim(n_dim);
  system_parameters params;
  params.n_dim = n_dim;
  params.system_radius = 50;
  params.boundary = BOX;
  params.graph_background = 1;
  params.delta = 0.0001;
  params.n_steps = 100000000;
  Object::SetParameters(&params);
  Object::SetDelta(params.delta);
  RNG rng(seed);

  // Init structures

  double pos[3] = {0,0,0};
  double u[3] = {1,0,0};
  Mesh m;
  m.Reserve(n_bonds);
  m.InitSiteAt(pos,1);

  for (int i=0; i<n_bonds; ++i) {
    double bl = 5+10*gsl_rng_uniform_pos(rng.r);
    //double bl = 10;
    m.AddRandomBondAnywhere(bl);
  }
  m.Report();
  fprintf(stderr,"----------------\n");
  m.SubReport();
  std::vector<Motor> mots;
  for (int i=0; i<n_motors; ++i) {
    Motor mot;
    mot.SetDiameter(2);
    mot.SetDiffusion();
    mots.push_back(mot);
    int i_bond = gsl_rng_uniform_int(rng.r, n_bonds);
    mots[i].AttachToBond(std::make_pair(m.GetBond(i_bond),NONE),m.GetBond(i_bond)->GetLength()*gsl_rng_uniform_pos(rng.r));
  }

   //Init space
  space_struct space;
  space.n_dim = params.n_dim;
  space.n_periodic = 0;
  space.radius = 3.0*params.system_radius;
  space.type = "box";
  double unit_cell[9] = {space.radius,0,0,0,space.radius,0,0,0,space.radius};
  space.unit_cell = unit_cell;

   //Init graphics

  double background_color = (params.graph_background == 0 ? 0.1 : 1);
  std::vector<graph_struct*> graph_array;
  m.Draw(&graph_array);
  for (int i=0; i<n_motors; ++i) {
    mots[i].Draw(&graph_array);
  }

  graphics.Init(&graph_array, &space, background_color, 0);
  graphics.DrawLoop();

   //Simulation loop

  int i_step=0;
  while (i_step<params.n_steps) {
    i_step++;

    // ...
    // Update simulation
    // ...

    for (int i=0; i<n_motors; ++i) {
      mots[i].UpdatePosition();
    }


    // update graphics
    if (i_step % 1000 == 0) {
      graph_array.clear();
      m.Draw(&graph_array);
      for (int i=0; i<n_motors; ++i) {
        mots[i].Draw(&graph_array);
      }
      graphics.Draw();
    }
  }

  // clean up (for when we actually terminate simulation)
  graphics.Clear();

  return 0;
}

