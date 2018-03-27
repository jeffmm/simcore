#ifndef _SIMCORE_ANALYSIS_MANAGER_H_
#define _SIMCORE_ANALYSIS_MANAGER_H_

#include "auxiliary.h"
#include "diffusion_analysis.h"
#include "filament_analysis.h"

class AnalysisManager {

  private:
    PositReader preader_;
    DiffusionAnalysis diff_analyzer_;
    FilamentAnalysis filament_analyzer_;

  public:
    void RunAnalyses(system_parameters params, std::vector<std::string> pfiles) {
      std::cout << "  Beginning analyses:\n";
      for (auto pfile : pfiles) {
        if (debug_trace) {
          preader_.PrintPosit(pfile);
        }
        if (params.diffusion_analysis) {
          std::cout << "    Performing diffusion analysis on " << pfile << "\n";
          diff_analyzer_.CalculateDiffusion(&params, pfile);
        }
        if (params.filament_analysis) {
          std::cout << "    Performing filament analysis on " << pfile << "\n";
          filament_analyzer_.AnalyzeFilament(&params, pfile);
        }
      }
      std::cout << "  Finished analyses.\n";
    }

};
#endif // _SIMCORE_ANALYSIS_MANAGER_H_
