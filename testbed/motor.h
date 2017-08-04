#include "bond.h"

class Motor : public Site {

  protected:
    bool bound_;
    double k_on_,
           k_off_,
           bond_length_,
           bond_lambda_,
           diffusion_;
  public:
    Motor(long seed);
    bool UpdatePriors();
    void UpdatePosition();
    void Diffuse();
    void DiffuseBound();
    void AttachToBond(directed_bond, double lambda);
    bool SwitchBonds(bool next_bond, double lambda);
    void UpdateMotorPosition() {}
    void SetDiffusion();
};

