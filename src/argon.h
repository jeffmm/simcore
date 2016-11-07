#ifndef _SIMCORE_ARGON_H_
#define _SIMCORE_ARGON_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

class Argon : public Simple {
  protected:
    double mass_;
  public:
    Argon(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {
      // Set parameters unique to MD bead
      diameter_ = params->argon_diameter;
      mass_ = params->argon_mass;
    }
    ~Argon() {}
    Argon(const Argon& that) : Simple(that) {}
    Argon& operator=(Argon const& that) {
      Simple::operator=(that); return *this;
    }
    virtual void Init();
    void InitConfigurator(std::array<double, 3> rx, std::array<double, 3> vx);
    virtual void UpdatePosition();
    virtual void UpdatePositionMP();
    virtual void Integrate();
    virtual void UpdateKineticEnergy();
    virtual double const GetKineticEnergy();

};

class ArgonSpecies : public Species<Argon> {
  protected:

  public:
    ArgonSpecies() : Species() {
      SetSID(SID::argon);
    }
    ArgonSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::argon);
    }
    ~ArgonSpecies() {}
    ArgonSpecies(const ArgonSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

    // Configurations
    void Configurator();
};

#endif // _SIMCORE_ARGON_H_
