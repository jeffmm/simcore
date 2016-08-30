#ifndef _SIMCORE_LOOKUP_TABLE_H_
#define _SIMCORE_LOOKUP_TABLE_H_

#include <string>
#include <vector>

#include "auxiliary.h"

class LookupTable {
  protected:
    int ndim_;

    std::vector<int> ngrid_;
    std::vector<int> ibin_;

    std::vector<double> *x_;
    std::vector<double> table_;

  public:
    void Init(int pndim, std::vector<double> *pX,
              double (*func)(std::vector<double> &x, void *params),
              void *params);
    double Lookup(double *x);
    double Lookup(double x, double y);
    double Invert(int dim, double u, double val[]);

    void OutputBinary(std::string outfile);
    void OutputInterpolatedGrid(std::vector<double> *x,
                                std::string outfile);

    // inlines
    inline int LinearIndex(int ibin) {
      return ibin_[0]*ngrid_[1] + ibin_[1];
    }

    inline void Lookup(std::vector<double> &pos) {
      Lookup(pos.data());
    }

};

#endif // _SIMCORE_LOOKUP_TABLE_H_
