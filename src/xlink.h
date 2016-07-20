#ifndef _SIMCORE_XLINK_H_
#define _SIMCORE_XLINK_H_

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
    double n_exp_;

    attach_type bound_;

    /* Functions unique to Br Dimer */
    void Integrate();
    void InternalForces();
    void UpdateOrientation();
    void CheckBoundState();
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
    void UpdatePositionMP();
    void Draw(std::vector<graph_struct*> * graph_array);
    void ApplyInteractions();
    void Dump() {
      Object::Dump();
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

    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    double const GetNExp() {return n_exp_;}
    attach_type const GetBoundState() {return bound_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
    void SetNExp(double const n) {n_exp_=n;}

    void DumpKMC() {
      auto head0 = elements_.begin();
      auto head1 = elements_.begin()+1;
      printf("\t\t[%d] -> {n_exp: %2.4f (%2.4f, %2.4f)}\n", GetOID(), n_exp_, head0->GetNExp(), head1->GetNExp());
    }
};

#include "species.h"
class XlinkSpecies : public Species<Xlink> {
  protected:
    double n_exp_ = 0.0;
    int nbound_ = 0.0;
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

    // KMC specifics
    double const GetNExp() {return n_exp_;}
    void SetNExp(double n) {n_exp_=n;}
    int const GetNFree() {return nfree_;}
    void SetNFree(int n) {nfree_=n;}
    int const GetNBound() {return nbound_;}
    void SetNBound(int n) {nbound_=n;}

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
