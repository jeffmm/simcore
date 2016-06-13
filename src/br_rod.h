#ifdef _SIMCORE_BR_ROD_H_
#define _SIMCORE_BR_ROD_H_

#include "br_bead.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

class BrRod : public Composite<BrBead,Simple> {

  private:
    double max_child_length_;
  public:
    BrRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Simple(params, space, seed, sid) {
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        max_child_length_ = 0.5*params->cell_length;
      }
    ~BrRod() {}
    BrRod(const BrRod& that) : Composite(that) {}
    BrRod& operator=(BrRod const& that) {Composite::operator=(that); return *this;} 
    virtual void Init() {Composie::Init()};
    virtual void UpdatePosition();
    virtual void Integrate();
};

class BrRodSpecies : public Species<BrRod> {
  protected:
    void InitPotentials(system_parameters *params);
  public:
    BrRodSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::br_rod);
      InitPotentials(params);
    }
    ~BrRodSpecies() {}
    BrRodSpecies(const BrRodSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

};

#endif // _SIMCORE_BR_ROD_H_
