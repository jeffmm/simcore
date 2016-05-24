#ifndef _CYTOSCORE_BROWNIAN_DIMER_H_
#define _CYTOSCORE_BROWNIAN_DIMER_H_

#include "auxiliary.h"
//#include "object.h"
#include "brownian_bead.h"

class BrownianDimer : public Composite<BrownianBead> {
  private:
    /* Unique to BrownianDimer */
    double eq_length_;
    double k_spring_;

    /* Functions unique to Brownian Dimer */
    void Integrate();
    void InternalForces();
    void UpdateOrientation();
    //void InsertRandom(); // Temporary, get space to do this
    graph_struct g2_; // Drawing two spheros, one for each bead;

  public:
    //Constructor
    BrownianDimer(system_parameters *params, space_struct *space, long seed, unsigned int const sid) : Composite(params, space, seed, sid) {
      for (int i=0; i<2; ++i) {
        BrownianBead b(params, space, gsl_rng_get(rng_.r), GetSID());
        b.SetCID(GetCID());
        elements_.push_back(b);
      }
      diameter_ = params->dimer_diameter;
      length_ = params->dimer_length;
      eq_length_ = params->dimer_eq_length;
      k_spring_ = params->dimer_k_spring;
    }
    //Destructor
    ~BrownianDimer() {}
    //Copy constructor
    BrownianDimer(const BrownianDimer& that) : Composite(that) {
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
    }
    //Assignment constructor
    BrownianDimer& operator=(BrownianDimer const& that) {
      Composite::operator=(that);
      eq_length_=that.eq_length_;
      k_spring_=that.k_spring_;
      g2_=that.g2_;
      return *this;
    };
    /* Define virtual functions */
    //void Init() {InsertRandom();} //TODO: Temporary, have space do this
    void Init();
    void UpdatePosition();
    void Draw(std::vector<graph_struct*> * graph_array);

    /* Functions unique to Brownian Dimer */
    double const GetKSpring() {return k_spring_;}
    double const GetEqLength() {return eq_length_;}
    void SetEqLength(double const eq_length) {eq_length_ = eq_length;}
    void SetKSpring(double const k_spring) {k_spring_ = k_spring;}
};

#include "species.h"
class BrownianDimerSpecies : public Species<BrownianDimer> {

  public:
    BrownianDimerSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {}
    ~BrownianDimerSpecies() {}
    BrownianDimerSpecies(const BrownianDimerSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
};

#endif // _CYTOSCORE_BROWNIAN_DIMER_H_
