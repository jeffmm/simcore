#ifndef _SIMCORE_XLINK_H_
#define _SIMCORE_XLINK_H_

#include "auxiliary.h"
#include "br_bead.h"

class Xlink : public Composite<BrBead> {
  private:
    /* Unique to Xlink */
    double eq_length_;
    double k_spring_;

    bool bound_;

    /* Functions unique to Br Dimer */
    void Integrate();
    void InternalForces();
    void UpdateOrientation();
    //void InsertRandom(); // Temporary, get space to do this
    graph_struct g2_; // Drawing two spheros, one for each bead;

  public:
    //Constructor
    Xlink(system_parameters *params, space_struct *space, long seed, SID sid) : Composite(params, space, seed, sid) {
      for (int i=0; i<2; ++i) {
        BrBead b(params, space, gsl_rng_get(rng_.r), GetSID());
        b.SetCID(GetCID());
        elements_.push_back(b);
        bound_ = false;
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
    void UpdatePosition();
    void UpdatePositionMP();
    void Draw(std::vector<graph_struct*> * graph_array);
    void ApplyInteractions();
    void Dump() {
      Object::Dump();
      printf("\tXlink {d: %2.2f}, {length: %2.2f}, {eq_length: %2.2f}, {k_spring: %2.2f}\n", diameter_, length_,
          eq_length_, k_spring_);
    }

    /* Functions unique to Xlink */
    void DiffuseXlink();

    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
};

#include "species.h"
class XlinkSpecies : public Species<Xlink> {
  protected:

  public:
    XlinkSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::xlink);
    }
    ~XlinkSpecies() {}
    XlinkSpecies(const XlinkSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _SIMCORE_XLINK_H_
