#ifndef _SIMCORE_XLINK_HEAD_H_
#define _SIMCORE_XLINK_HEAD_H_

#include "object.h"
#include "auxiliary.h"

#include <iomanip>

class XlinkHead : public Simple {
  private:
    double diffusion_;

    // kmc stuff
    double n_exp_0_1_ = 0.0;
    double n_exp_1_2_ = 0.0;
    bool bound_;
    int attachidx_ = -1;
    double attachpos_ = 0.0;
  public:
    XlinkHead(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {
      diameter_=params->br_walker_diameter;
      is_kmc_ = true;
      bound_ = false;
      SetDiffusion();
    }
    ~XlinkHead() {}
    XlinkHead(const XlinkHead& that) : Simple(that) {}
    XlinkHead& operator=(XlinkHead const& that) {Simple::operator=(that); return *this;} 
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePositionMP();
    void Init();

    // kmc specifics
    virtual void PrepKMC(std::vector<neighbor_t>* neighbors);
    virtual void StepKMC();
    double const GetNExp_0_1() {return n_exp_0_1_;}
    double const GetNExp_1_2() {return n_exp_1_2_;}
    void SetNExp_0_1(double const n) {n_exp_0_1_=n;}
    void SetNExp_1_2(double const n) {n_exp_1_2_=n;}
    bool const GetBound() {return bound_;}
    void SetBound(bool bound) {bound_=bound;}
    void Attach(int idx, double pos);
    std::pair<int, double> GetAttach() {return std::make_pair(attachidx_, attachpos_);}

    virtual void DumpKMC() {
      if (n_exp_0_1_ < 0.0) {
        printf("ERROR - nexp_0_1 < 0.0\n");
        exit(1);
      }
      if (n_exp_1_2_ < 0.0) {
        printf("ERROR - nexp_1_2 < 0.0\n");
        exit(1);
      }
      std::cout << std::setprecision(16) << "          head[" << GetOID() << "] -> {n_exp: ";
      std::cout << std::setprecision(16) << "(" << GetNExp_0_1() << ", " << GetNExp_1_2();
      std::cout << std::setprecision(16) << ")}, {bound: " << (GetBound() ? "true " : "false") << "}, ";
      std::cout << std::setprecision(16) << "{pidx: " << GetAttach().first << ", " << GetAttach().second;
      std::cout << std::setprecision(16) << "}\n";
    }
};

#include "species.h"
class XlinkHeadSpecies : public Species<XlinkHead> {
  protected:
    double n_exp_0_1_ = 0.0;
    double n_exp_1_2_ = 0.0;
    int nbound_ = 0.0;
    int nfree_ = 0.0;

  public:
    XlinkHeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::xlink_head);
      is_kmc_ = true;
    }
    ~XlinkHeadSpecies() {}
    XlinkHeadSpecies(const XlinkHeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    virtual void DumpKMC() {
      for (auto it = members_.begin(); it != members_.end(); ++it) {
        printf("\t\t[%d] -> {n_exp: (%2.4f, %2.4f)}, {bound: %s}, {parent index: %d, %2.4f}\n", (*it)->GetOID(), (*it)->GetNExp_0_1(),
            (*it)->GetNExp_1_2(), (*it)->GetBound() ? "true " : "false",
            (*it)->GetAttach().first, (*it)->GetAttach().second);
      }
    }


    double const GetNExp_0_1() {return n_exp_0_1_;}
    double const GetNExp_1_2() {return n_exp_1_2_;}
    void SetNExp_0_1(double n) {n_exp_0_1_=n;}
    void SetNExp_1_2(double n) {n_exp_1_2_=n;}
    int const GetNFree() {return nfree_;}
    void SetNFree(int n) {nfree_=n;}
    int const GetNBound() {return nbound_;}
    void SetNBound(int n) {nbound_=n;}
};

#endif // _SIMCORE_XLINK_HEAD_H_

