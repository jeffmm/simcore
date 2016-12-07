#ifndef _SIMCORE_BR_SIMPLE_ROD_H_
#define _SIMCORE_BR_SIMPLE_ROD_H_

#include "object.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

class BrSimpleRod : public Simple {
  private:
    double gamma_par_,
           gamma_perp_,
           gamma_rot_,
           rand_sigma_par_,
           rand_sigma_perp_,
           rand_sigma_rot_,
           body_frame_[6];
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
  public:
    BrSimpleRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Simple(params, space, seed, sid) {
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        double max_length = 0.5*params->cell_length;
        if (length_ > max_length)
          length_=max_length;
      }
    void Init();
    void UpdatePosition();
    void UpdatePositionMP();
    void Integrate();
};

class BrSimpleRodSpecies : public Species<BrSimpleRod> {
  protected:
    //void InitPotentials(system_parameters *params);
  public:
    BrSimpleRodSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::br_simple_rod);
      //InitPotentials(params);
    }
    void Init() {
      Species::Init();
    }

};

#endif // _SIMCORE_BR_SIMPLE_ROD_H_
