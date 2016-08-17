#ifndef _SIMCORE_XLINK_H_
#define _SIMCORE_XLINK_H_

#include <cassert>

#include "auxiliary.h"
#include "xlink_head.h"

enum attach_type {
  unbound = 0,
  singly,
  doubly
};

class Xlink : public Composite<XlinkHead> {
  private:
    /* Unique to Xlink */
    double eq_length_;
    double k_spring_;
    double n_exp_0_1_ = 0.0;
    double n_exp_1_2_ = 0.0;
    double uinternal_ = 0.0;

    attach_type bound_;

    /* Functions unique to Br Dimer */
    void Integrate();
    void UpdateOrientation();
    graph_struct g2_; // Drawing two spheros, one for each bead;

  public:
    //Constructor
    Xlink(system_parameters *params, space_struct *space, long seed, SID sid) : Composite(params, space, seed, sid) {
      for (int i=0; i<2; ++i) {
        XlinkHead b(params, space, gsl_rng_get(rng_.r), SID::xlink_head);
        b.SetCID(GetCID());
        elements_.push_back(b);
        bound_ = unbound;
        is_kmc_ = true;
      }
      diameter_ = params->xlink_diameter;
      length_ = 0.0;
      eq_length_ = params->dimer_eq_length;
      k_spring_ = params->dimer_k_spring;
    }
    //Destructor
    ~Xlink() {}
    //Copy constructor
    Xlink(const Xlink& that) : Composite(that) {
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
    }
    //Assignment constructor
    Xlink& operator=(Xlink const& that) {
      Composite::operator=(that);
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
      return *this;
    };
    /* Define virtual functions */
    void Init();
    void InitConfigurator(const double* const x, const double diameter);
    void UpdatePositionMP();
    void Draw(std::vector<graph_struct*> * graph_array);
    void ApplyInteractions();
    void Dump() {
      Object::Dump();
      if (isnan(position_[0])) {
        printf("Something has gone horribly wrong!\n");
        exit(1);
      }
      printf("\t   {d: %2.2f}, {length: %2.2f}, {eq_length: %2.2f}, {k_spring: %2.2f}, {bound: %d}\n", diameter_, length_,
          eq_length_, k_spring_, (int)bound_);
      auto head0 = elements_.begin();
      auto head1 = elements_.begin()+1;
      printf("\t   ");
      head0->Dump();
      printf("\t   ");
      head1->Dump();
    }

    /* Functions unique to Xlink */
    void DiffuseXlink();
    std::vector<XlinkHead>* GetHeads() {return &elements_;}
    std::pair<bool,bool> GetBoundHeads(XlinkHead **freehead, XlinkHead **boundhead);
    void CheckBoundState();

    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    double const GetNExp_0_1() {return n_exp_0_1_;}
    double const GetNExp_1_2() {return n_exp_1_2_;}
    double const GetInternalU() {return uinternal_;}
    attach_type const GetBoundState() {return bound_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
    void SetNExp_0_1(double const n) {n_exp_0_1_=n;}
    void SetNExp_1_2(double const n) {n_exp_1_2_=n;}
    void SetInternalU(double u) {uinternal_=u;}

    void DumpKMC() {
      auto head0 = elements_.begin();
      auto head1 = elements_.begin()+1;
      assert(n_exp_0_1_ < 10000);
      assert(n_exp_1_2_ < 10000);
      printf("\t\t[%d] -> {u: %2.4f}, {n_exp_0_1: %2.4f (%2.4f, %2.4f)}, {n_exp_1_2: %2.4f (%2.4f, %2.4f)}\n", GetOID(),
          uinternal_,
          n_exp_0_1_, head0->GetNExp_0_1(), head1->GetNExp_0_1(),
          n_exp_1_2_, head0->GetNExp_1_2(), head1->GetNExp_1_2());
      head0->DumpKMC();
      head1->DumpKMC();
    }
};

#include "species.h"
class XlinkSpecies : public Species<Xlink> {
  protected:
    double n_exp_0_1_ = 0.0;
    double n_exp_1_2_ = 0.0;
    int nbound1_[2] = {0, 0};
    int nbound2_[2] = {0, 0};
    int nfree_ = 0.0;

  public:
    XlinkSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::xlink);
      is_kmc_ = true;
    }
    ~XlinkSpecies() {}
    XlinkSpecies(const XlinkSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }

    // Specific initialization
    void ConfiguratorXlink();

    // KMC specifics
    double const GetNExp_0_1() {return n_exp_0_1_;}
    double const GetNExp_1_2() {return n_exp_1_2_;}
    void SetNExp_0_1(double n) {n_exp_0_1_=n;}
    void SetNExp_1_2(double n) {n_exp_1_2_=n;}
    int const GetNFree() {return nfree_;}
    void SetNFree(int n) {nfree_=n;}
    const int* const GetNBound1() {return nbound1_;}
    void SetNBound1(int n0, int n1) {nbound1_[0]=n0; nbound1_[1]=n1;}
    const int* const GetNBound2() {return nbound2_;}
    void SetNBound2(int n0, int n1) {nbound2_[0]=n0; nbound2_[1]=n1;}

    virtual void DumpKMC() {
      for (auto it = members_.begin(); it != members_.end(); ++it) {
        (*it)->DumpKMC();
      }
    }

    std::vector<Xlink*>* GetXlinks() {
      return &members_;
    }
};

#endif // _SIMCORE_XLINK_H_
