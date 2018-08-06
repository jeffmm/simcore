#ifndef _SIMCORE_STRUCT_ANALYSIS_H_
#define _SIMCORE_STRUCT_ANALYSIS_H_

#include "auxiliary.h"
#include "interaction.h"
#include "object.h"

class StructAnalysis {
  private:
    system_parameters * params_;
    int count_,
        n_dim_,
        n_bins_,
        n_bins_1d_,
        average_structure_,
        n_objs_;
    double bin_width_;
    float * pdf_array_,
          * nematic_array_,
          * polar_array_,
          * pdf_array_temp_,
          * nematic_array_temp_,
          * polar_array_temp_;
    std::fstream pdf_file_;
    std::fstream nematic_file_;
    std::fstream polar_file_;

  public:
    void Init(system_parameters *params);
    void Clear();
    void CalculateStructurePair(std::vector<pair_interaction>::iterator pix);
    void WriteStructData();
    void PopulateLine(double * site1, double * site2, double dotprod);
    void BinLine(int x0, int y0, int x1, int y1, double dotprod);
    void BinLineLow(int x0, int y0, int x1, int y1, double dotprod);
    void BinLineHigh(int x0, int y0, int x1, int y1, double dotprod);
    void BinArray(int x, int y, double dotprod);
    void AverageStructure();
    void IncrementCount() {count_++;}
    void SetNumObjs(int nobj) {n_objs_ = nobj;}
    int GetNumObjs() {return n_objs_;}
};

#endif
