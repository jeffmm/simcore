#include "object.h"
#include "omp.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>

class Simulation {
  private:
    std::vector<Object> objs;
    void InitRNG();
    void InitPositions();
    void CreatePairs();
  public:
    void Run();
};
