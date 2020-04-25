#include "simcore/potential_base.hpp"

double PotentialBase::max_force_ = -1;
int PotentialBase::n_dim_ = -1;
bool PotentialBase::fmax_violation_ = false;
void PotentialBase::SetMaxForce(double fmax) { max_force_ = fmax; }
double PotentialBase::GetMaxForce() { return max_force_; }
bool PotentialBase::CheckMaxForceViolation() { return fmax_violation_; }
void PotentialBase::ResetMaxForceViolation() { fmax_violation_ = false; }
void PotentialBase::MaxForceViolation() { fmax_violation_ = true; }
void PotentialBase::SetNDim(int ndim) { n_dim_ = ndim; }
void PotentialBase::Init(system_parameters *params) {
  SetMaxForce(params->f_cutoff);
  SetNDim(params->n_dim);
  InitPotentialParams(params);
}
