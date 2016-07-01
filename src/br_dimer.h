#ifndef _SIMCORE_BR_DIMER_H_
#define _SIMCORE_BR_DIMER_H_

#include "auxiliary.h"
//#include "object.h"
#include "br_bead.h"
#include "lennard_jones_12_6.h"

class BrDimer : public Composite<BrBead> {
  private:
    /* Unique to BrDimer */
    double eq_length_;
    double k_spring_;

    /* Functions unique to Br Dimer */
    void Integrate();
    void InternalForces();
    void UpdateOrientation();
    //void InsertRandom(); // Temporary, get space to do this
    graph_struct g2_; // Drawing two spheros, one for each bead;

  public:
    //Constructor
    BrDimer(system_parameters *params, space_struct *space, long seed, SID sid) : Composite(params, space, seed, sid) {
      for (int i=0; i<2; ++i) {
        BrBead b(params, space, gsl_rng_get(rng_.r), GetSID());
        b.SetCID(GetCID());
        elements_.push_back(b);
      }
      diameter_ = params->dimer_diameter;
      length_ = params->dimer_length;
      eq_length_ = params->dimer_eq_length;
      k_spring_ = params->dimer_k_spring;
    }
    //Destructor
    ~BrDimer() {}
    //Copy constructor
    BrDimer(const BrDimer& that) : Composite(that) {
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
    }
    //Assignment constructor
    BrDimer& operator=(BrDimer const& that) {
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

    /* Functions unique to Br Dimer */
    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
};

#include "species.h"
class BrDimerSpecies : public Species<BrDimer> {
  protected:
    //void InitPotentials (system_parameters *params);

  public:
    BrDimerSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_dimer);
      //InitPotentials(params);
    }
    ~BrDimerSpecies() {}
    BrDimerSpecies(const BrDimerSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _SIMCORE_BR_DIMER_H_
