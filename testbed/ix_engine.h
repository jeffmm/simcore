#include "particle_tracker.h"

class IxEngine {
  private:
    ParticleTracker pt;
    std::vector<Object> * objs;
    std::vector<ix_pair> nlist;
    void Interact(int i, int j, int *count);
    bool CheckUpdate();
    void UpdatePos0s();
  public:
    IxEngine(std::vector<Object> * o);
    void Init();
    void TrackParticles();
    void InteractPairs();
    void Cleanup();
};
