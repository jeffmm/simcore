#ifndef _SIMCORE_ANALYSIS_MANAGER_H_
#define _SIMCORE_ANALYSIS_MANAGER_H_

#include "auxiliary.h"
#include "diffusion_analysis.h"

class AnalysisManager {

  private:
    DiffusionAnalysis diff_analyzer_;

  public:
    void RunAnalyses(system_parameters params, std::vector<std::string> pfiles) {
      for (auto pfile : pfiles) {
        if (params.diffusion_analysis)
          diff_analyzer_.CalculateDiffusion(pfile);
      }
    }

};
#endif // _SIMCORE_ANALYSIS_MANAGER_H_
