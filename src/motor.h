#include "bond.h"

class Motor : public Site {

  protected:
    bool bound_,
         walker_;
    int step_direction_;
    double k_on_,
           k_off_,
           bond_length_,
           bond_lambda_,
           diffusion_;
  public:
    Motor();
    void Init();
    bool UpdatePriors();
    void UpdatePosition();
    void Diffuse();
    void DiffuseBound();
    void AttachToBond(directed_bond, double lambda);
    bool SwitchBonds(bool next_bond, double lambda);
    void UpdateMotorPosition() {}
    void SetDiffusion();
    void SetWalker(int dir);
};

