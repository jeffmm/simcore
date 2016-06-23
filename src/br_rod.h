#ifndef _SIMCORE_BR_ROD_H_
#define _SIMCORE_BR_ROD_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

class BrRod : public Composite<Site,Bond> {

  private:
    int n_bonds_;
    double max_child_length_,
           child_length_,
           gamma_par_,
           gamma_perp_,
           gamma_rot_,
           rand_sigma_par_,
           rand_sigma_perp_,
           rand_sigma_rot_,
           body_frame_[6];
    void UpdateSiteBondPositions();
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void ApplyForcesTorques();

  public:
    BrRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Composite(params, space, seed, sid) {
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        max_child_length_ = 0.5*params->cell_length;
        // Initialize end sites
        for (int i=0; i<2; ++i) {
          Site s(params, space, gsl_rng_get(rng_.r), GetSID());
          s.SetCID(GetCID());
          elements_.push_back(s);
        }
        // Initialize bonds
        n_bonds_ = (int) ceil(length_/max_child_length_);
        child_length_ = length_/n_bonds_;
        for (int i=0; i<n_bonds_; ++i) {
          Bond b(params, space, gsl_rng_get(rng_.r), GetSID());
          b.SetCID(GetCID());
          v_elements_.push_back(b);
        }
      }
    ~BrRod() {}
    BrRod(const BrRod& that) : Composite(that) {}
    BrRod& operator=(BrRod const& that) {Composite::operator=(that); return *this;} 
    virtual void Init();
    virtual void Integrate();
    virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositionMP();
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
