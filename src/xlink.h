#ifndef _SIMCORE_XLINK_H_
#define _SIMCORE_XLINK_H_

#include <cassert>

#include "auxiliary.h"
#include "xlink_head.h"

#include <iomanip>
#include <unordered_map>

enum attach_type {
  unbound = 0,
  singly,
  doubly
};

class Xlink : public Composite<XlinkHead> {
  private:
    /* Unique to Xlink */
    double n_exp_0_1_ = 0.0;
    double n_exp_1_2_ = 0.0;
    double r_cross_[3] = {0.0, 0.0, 0.0};

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
    }
    //Destructor
    ~Xlink() {}
    //Copy constructor
    Xlink(const Xlink& that) : Composite(that) {
      g2_=that.g2_;
    }
    //Assignment constructor
    Xlink& operator=(Xlink const& that) {
      Composite::operator=(that);
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
      std::cout << std::setprecision(16) << "{" << GetOID() << "," << GetRID() << "," << GetCID() << "}"
        << " -> x(" << GetPosition()[0] << ", " << GetPosition()[1] << ", " << GetPosition()[2] << "), "
        << "f(" << GetForce()[0] << ", " << GetForce()[1] << ", " << GetForce()[2] << "), "
        << "u(" << GetKineticEnergy() << "), p(" << GetPotentialEnergy() << "), "
        << "d(" << diameter_ << "), l(" << length_ << ")\n";
      auto head0 = elements_.begin();
      auto head1 = elements_.begin()+1;
      std::cout << "\theads ->\n";
      std::cout << "\t  ";
      head0->Dump();
      std::cout << "\t  ";
      head1->Dump();
    }
    std::vector<std::pair<unsigned int, unsigned int>> GetInternalPairs() {
      auto head0 = elements_.begin();
      auto head1 = elements_.begin() + 1;
      std::pair<unsigned int, unsigned int> mpair= std::make_pair(head0->GetOID(), head1->GetOID());
      std::vector<std::pair<unsigned int, unsigned int>> retval;
      retval.push_back(mpair);
      return retval;
    }

    /* Functions unique to Xlink */
    void DiffuseXlink();
    std::vector<XlinkHead>* GetHeads() {return &elements_;}
    std::pair<bool,bool> GetBoundHeads(XlinkHead **freehead, XlinkHead **boundhead);
    void CheckBoundState();

    void BindHeadSingle(int ihead, double crosspos, int rodoid);
    void BindHeadDouble(double crosspos0, int rodoid0, double crosspos1, int rodoid1);

    void UpdateStagePosition(const double* const xr0, const double* const ur0, const double lr0, const int atidx0,
                             const double* const xr1, const double* const ur1, const double lr1, const int atidx1);
    void UpdateStage1Position(const double* const xr0, const double* const ur0, const double lr0, const int atidx);
    void UpdateStage2Position(const double* const xr0, const double* const ur0, const double lr0, const int atidx0,
                             const double* const xr1, const double* const ur1, const double lr1, const int atidx1);

    double const GetNExp_0_1() {return n_exp_0_1_;}
    double const GetNExp_1_2() {return n_exp_1_2_;}
    attach_type const GetBoundState() {return bound_;}
    void SetNExp_0_1(double const n) {n_exp_0_1_=n;}
    void SetNExp_1_2(double const n) {n_exp_1_2_=n;}
    const double GetInternalEnergy();

    const double* const GetRcross() {return r_cross_;}

    void DumpKMC() {
      auto head0 = elements_.begin();
      auto head1 = elements_.begin()+1;
      assert(n_exp_0_1_ < 10000);
      assert(n_exp_1_2_ < 10000);

      std::string boundstring;
      switch(bound_) {
        case unbound:
          boundstring = "unbound";
          break;
        case singly:
          boundstring = "singly";
          break;
        case doubly:
          boundstring = "doubly";
          break;
      }
      std::cout << std::setprecision(16) << "        [" << GetOID() << "] -> {" << boundstring << "}";
      std::cout << std::setprecision(16) << ", {n_exp_0_1: " << n_exp_0_1_ << " (" << head0->GetNExp_0_1() << ", ";
      std::cout << std::setprecision(16) << head1->GetNExp_0_1() << ")}, {n_exp_1_2: " << n_exp_1_2_ << " (";
      std::cout << std::setprecision(16) << head0->GetNExp_1_2() << ", " << head1->GetNExp_1_2() << ")}\n";
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
    XlinkSpecies() : Species() {
      SetSID(SID::xlink);
      is_kmc_ = true;
    }
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
    void Configurator();

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

    static void CreateTestXlink(Xlink **mxit,
                                int ndim,
                                std::vector<Simple*>* simples,
                                std::unordered_map<int, int>* oid_position_map,
                                const std::string &filename,
                                const std::string &modulename,
                                const std::string &unitname,
                                const std::string &xname,
                                int itest,
                                int attachoid);

    static void CreateTestXlink(Xlink **mxit,
                                int ndim,
                                std::vector<Simple*>* simples,
                                std::unordered_map<int, int>* oid_position_map,
                                const std::string &filename,
                                const std::string &modulename,
                                const std::string &unitname,
                                const std::string &xname,
                                int itest,
                                int attachoid0,
                                int attachoid1);
};

#endif // _SIMCORE_XLINK_H_
