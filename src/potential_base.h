#ifndef _CYTOSCORE_POTENTIAL_BASE_H_
#define _CYTOSCORE_POTENTIAL_BASE_H_

// Potential base class for all external(or maybe even internal) potentials
class PotentialBase {
  protected:
    int ndim_;
    double rcut_, rcut2_; // Cutoff radius
    double box_, boxby2_; // Unit cell size (and half!)
  public:
    PotentialBase(int pNdim, double pRcut, double pBox) : ndim_(pNdim), rcut_(pRcut), box_(pBox) {
        rcut2_ = rcut_ * rcut_;
        boxby2_ = 0.5 * box_;
    }
    virtual ~PotentialBase() {}
    virtual void CalcPotential(double* x,
                               double* y,
                               double* fpote) {}
    virtual void Print() {
        std::cout << "\t{rcut:" << rcut_ << "}, {rcut2_:" << rcut2_ << "}, {box:" << box_ << "}\n";
    }
};

#endif
