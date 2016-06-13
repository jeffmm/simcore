#ifndef _SIMCORE_POTENTIAL_BASE_H_
#define _SIMCORE_POTENTIAL_BASE_H_

#include "auxiliary.h"

// Potential base class for all external(or maybe even internal) potentials
class PotentialBase {
  protected:
    int n_dim_;
    double rcut_, rcut2_; // Cutoff radius
    std::string pot_name_;
    space_struct *space_;
  public:
    PotentialBase(space_struct* pSpace, double pRcut) : rcut_(pRcut), space_(pSpace) {
      pot_name_ = "Untitled";
      rcut2_ = rcut_ * rcut_;
      n_dim_ = space_->n_dim;
    }
    virtual ~PotentialBase() {}
    virtual void CalcPotential(double *dr,
                       double dr_mag,
                       double buffer,
                       double *fpote) {}
    double GetRCut() {return rcut_;}
    double GetRCut2() {return rcut2_;}
    virtual void Print() {
      std::cout << pot_name_ << "\n";
        printf("\t{rcut:%2.2f}\n", rcut_);
        //for (int i = 0; i < n_dim_; ++i) {
            //printf("\t");
            //for (int j = 0; j < n_dim_; ++j) {
                //printf("%2.2f\t", space_->unit_cell[i][j]);
            //}
            //printf("\n");
        //}
    }
};

typedef std::pair<sid_pair, PotentialBase*> potential_pair;
typedef std::map<sid_pair, PotentialBase*> potential_map;

#endif // _SIMCORE_POTENTIAL_BASE_H_
