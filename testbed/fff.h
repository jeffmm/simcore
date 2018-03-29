#include "site.h"

class FFF {
  private:
    site_it filament_begin;
    site_it filament_end;
    void ApplyForcesTorques();
    void Integrate();
  public:
    void UpdatePosition();
};

class FFFEnsemble {
  private:
    midstep_;
    std::vector<Site> sites_;
    std::vector<FFF> fffs_;
    void AllocateFilaments();
  public:
    void UpdatePositions();
};

