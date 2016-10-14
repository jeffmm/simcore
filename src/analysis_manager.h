#ifndef _SIMCORE_ANALYSIS_MANAGER_H_
#define _SIMCORE_ANALYSIS_MANAGER_H_

#include "auxiliary.h"
#include "diffusion_analysis.h"

class AnalysisManager {

  private:
    DiffusionAnalysis diff_analyzer_;

  public:
    void RunAnalyses(system_parameters params, std::vector<std::string> pfiles) {
      std::cout << "  Beginning analyses:\n";
      for (auto pfile : pfiles) {
        if (params.diffusion_analysis) {
          std::cout << "    Performing diffusion analysis on " << pfile << "\n";
          diff_analyzer_.CalculateDiffusion(&params, pfile);
        }
      }
      std::cout << "  Finished analyses.\n";
    }

};
#endif // _SIMCORE_ANALYSIS_MANAGER_H_
