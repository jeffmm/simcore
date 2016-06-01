#ifndef _CYTOSCORE_TEST_POTENTIAL_H_
#define _CYTOSCORE_TEST_POTENTIAL_H_

#include "auxiliary.h"
#include "potential_base.h"

// Potential base class for all external(or maybe even internal) potentials
class TestPotential : public PotentialBase {
  public:
    TestPotential(space_struct *pSpace, double pRcut) : PotentialBase(pSpace, pRcut) {
      pot_name_ = "Test Potential";
    }
    virtual ~TestPotential() {}
    void CalcPotential(double* dr,
                       double dr_mag,
                       double buffer,
                       double* fpot) {
      double mag2_inv = dr_mag-buffer;
      mag2_inv = 1.0/(mag2_inv*mag2_inv);
      for (int i=0; i<n_dim_; ++i) {
        fpot[i] = - mag2_inv * dr[i] / dr_mag;
      }
    }

};

//class Interaction {
  //private:
    //SID sid1_,
        //sid2_;
    //PotentialBase * potential_;
  //public:
    //Interaction(SID sid1, SID sid2, PotentialBase *pot) :
      //sid1_(sid1), sid2_(sid2), potential_(pot) {}
    //PotentialBase * GetPotential() {return potential_;}
    //SID const GetSID1() {return sid1_;}
    //SID const GetSID2() {return sid2_;}
//};

#endif
