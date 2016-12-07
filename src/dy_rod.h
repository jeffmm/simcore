#ifndef _SIMCORE_DY_ROD_H_
#define _SIMCORE_DY_ROD_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

#include <unordered_map>

class DyRod : public Composite<Site,Bond> {

  private:
    int n_bonds_,
        diffusion_validation_flag_;
    double max_length_,
           min_length_,
           max_child_length_,
           child_length_,
           gamma_par_,
           gamma_perp_,
           gamma_rot_,
           rand_sigma_par_,
           rand_sigma_perp_,
           rand_sigma_rot_,
           body_frame_[6],
           driving_factor_;
    poly_state_t poly_state_;
    void UpdateSitePositions();
    void UpdateBondPositions();
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void ApplyForcesTorques();

  public:
    DyRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Composite(params, space, seed, sid) {
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        max_length_ = params->max_rod_length;
        min_length_ = params->min_rod_length;
        max_child_length_ = 0.5*params->cell_length;
        driving_factor_ = params->driving_factor;
        diffusion_validation_flag_ = params->diffusion_validation_flag;
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
          b.SetRID(GetRID());
          v_elements_.push_back(b);
        }
      }
    virtual void Init();
    virtual void Integrate();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition();
    virtual void Dump();

    // Specific functions for configurations
    void InitConfigurator(const double* const x, const double* const u, const double l);
};

class DyRodSpecies : public Species<DyRod> {
  protected:
    //void InitPotentials(system_parameters *params);
    double max_length_;
    double min_length_;
  public:
    DyRodSpecies() : Species() {
      SetSID(SID::dy_rod);
    }
    DyRodSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::dy_rod);
      //InitPotentials(params);
      max_length_ = params_->max_rod_length;
    }
    void Init() {
      Species::Init();
    }
    double const GetMaxLength() {return max_length_;}
    double const GetMinLength() {return min_length_;}

    void UpdatePositions() {
      for (auto it=members_.begin(); it != members_.end(); ++it) {
        (*it)->UpdatePosition();
      }
    }
    // Special insertion routine
    void Configurator();
};

#endif // _SIMCORE_DY_ROD_H_
