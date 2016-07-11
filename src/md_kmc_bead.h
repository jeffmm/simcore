#ifndef _SIMCORE_MD_KMC_BEAD_H_
#define _SIMCORE_MD_KMC_BEAD_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"
#include "lennard_jones_12_6.h"
#include "wca.h"

class MDKMCBead : public Simple {
  protected:
    double mass_;

    // kmc stuff
    double n_exp_; // number of expected for this particular bead this timestep
    double eps_eff_ = 30910;
    double on_rate_ = 0.003916;
    bool bound_;
    int attachidx_ = -1;
  public:
    MDKMCBead(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {
      // Set parameters unique to MD KMC bead
      diameter_ = params->md_kmc_bead_diameter;
      mass_ = params->md_kmc_bead_mass;
      is_kmc_ = true;
      bound_ = false;
    }
    ~MDKMCBead() {}
    MDKMCBead(const MDKMCBead& that) : Simple(that) {}
    MDKMCBead& operator=(MDKMCBead const& that) {
      Simple::operator=(that); return *this;
    }
    virtual void Init();
    virtual void UpdatePosition();
    virtual void UpdatePositionMP();
    virtual void Integrate();
    virtual void UpdateKineticEnergy();
    virtual double const GetKineticEnergy();

    // kmc specifics
    virtual void PrepKMC(std::vector<neighbor_t>* neighbors);
    virtual void StepKMC();
    double const GetNExp() {return n_exp_;}
    void SetNExp(double n) {n_exp_=n;}
    bool const GetBound() {return bound_;}
    void SetBound(bool bound) {bound_=bound;}
    void Attach(int idx);
    int const GetAttach() {return attachidx_;}
    void SelfSetPosition();
};

class MDKMCBeadSpecies : public Species<MDKMCBead> {
  protected:
    //void InitPotentials (system_parameters *params);
    double n_exp_;
    int nbound_;
    int nfree_;
  public:
    MDKMCBeadSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::md_kmc_bead);
      //InitPotentials(params);
      is_kmc_ = true;
      nbound_ = 0;
      n_exp_ = 0.0;
    }
    ~MDKMCBeadSpecies() {}
    MDKMCBeadSpecies(const MDKMCBeadSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

    virtual void PrepKMC();
    virtual void StepKMC();
    virtual void DumpKMC();
    int& GetNBound() {return nbound_;}
    int& GetNFree() {return nfree_;}
};

#endif // _SIMCORE_MD_KMC_BEAD_H_
