#include "ix_engine.h"

class FFFSim {
  private:
    IxEngine ix_engine;
    std::vector<Object> objs;
    void InitRNG();
    void InitPositions();
    void UpdatePositions();
    //void CreatePairs();
  public:
    FFFSim() : ix_engine(&objs) {}
    void Run();
};
