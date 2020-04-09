#include "simcore/analysis.hpp"

void AnalysisBase::SetParams(system_parameters *params) { params_ = params; }
void AnalysisBase::SetSpace(space_struct *space) { space_ = space; }
const space_struct *AnalysisBase::space_ = nullptr;
const system_parameters *AnalysisBase::params_ = nullptr;

