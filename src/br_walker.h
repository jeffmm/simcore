#ifndef _SIMCORE_BR_WALKER_H_
#define _SIMCORE_BR_WALKER_H_

#include "object.h"
#include "auxiliary.h"
#include "wca.h"

class BrWalker : public Simple {
  private:
    double diffusion_;

    // kmc stuff
    double n_exp_;
    bool bound_;
    int attachidx_ = -1;
    double attachpos_ = 0.0;
  public:
    BrWalker(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {
      diameter_=params->br_walker_diameter;
      is_kmc_ = true;
      bound_ = false;
      SetDiffusion();
    }
    ~BrWalker() {}
    BrWalker(const BrWalker& that) : Simple(that) {}
    BrWalker& operator=(BrWalker const& that) {Simple::operator=(that); return *this;} 
    void SetDiffusion() {diffusion_ = sqrt(24.0*diameter_/delta_);}
    void KickBead();
    void UpdatePosition();
    void UpdatePositionMP();
    void Init();

    // kmc specifics
    virtual void PrepKMC(std::vector<neighbor_t>* neighbors);
    virtual void StepKMC();
    double const GetNExp() {return n_exp_;}
    void SetNExp(double n) {n_exp_=n;}
    bool const GetBound() {return bound_;}
    void SetBound(bool bound) {bound_=bound;}
    void Attach(int idx, double pos);
    std::pair<int, double> GetAttach() {return std::make_pair(attachidx_, attachpos_);}
};

#include "species.h"
class BrWalkerSpecies : public Species<BrWalker> {
  protected:
    //void InitPotentials(system_parameters *params);
    double n_exp_ = 0.0;
    int nbound_ = 0.0;
    int nfree_ = 0.0;

  public:
    BrWalkerSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_walker);
      is_kmc_ = true;
      //InitPotentials(params);
    }
    ~BrWalkerSpecies() {}
    BrWalkerSpecies(const BrWalkerSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    virtual void DumpKMC() {
      for (auto it = members_.begin(); it != members_.end(); ++it) {
        printf("\t\t[%d] -> {n_exp: %2.4f}, {bound: %s}, {parent index: %d, %2.4f}\n", (*it)->GetOID(), (*it)->GetNExp(), (*it)->GetBound() ? "true " : "false",
            (*it)->GetAttach().first, (*it)->GetAttach().second);
      }
    }


    double const GetNExp() {return n_exp_;}
    void SetNExp(double n) {n_exp_=n;}
    int const GetNFree() {return nfree_;}
    void SetNFree(int n) {nfree_=n;}
    int const GetNBound() {return nbound_;}
    void SetNBound(int n) {nbound_=n;}
};

#endif // _SIMCORE_BR_WALKER_H_

