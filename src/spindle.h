#ifndef _SIMCORE_SPINDLE_H_
#define _SIMCORE_SPINDLE_H_

#include "spherocylinder.h"
#include "filament.h"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class Spindle : public Spherocylinder {
  protected:
    bool alignment_potential_,
         fixed_spacing_;
    int n_filaments_bud_,
        n_filaments_mother_;
    double k_spring_,
           k_align_,
           spring_length_,
           anchor_distance_,
           gamma_trans_,
           gamma_rot_,
           diffusion_,
           spb_diameter_;
    std::vector<Filament> filaments_;
    std::vector<Anchor> anchors_;
    void ApplyForcesTorques();
    void ApplyBoundaryForces();
    void InsertSpindle();
    void GenerateAnchorSites();
    void Integrate();
    void ResetAnchorPositions();
    bool InsertFilament(int i);
  public:
    Spindle();
    void Init();
    void UpdatePosition() {}
    void UpdatePosition(bool midstep);
    virtual std::vector<Object*> GetInteractors();
    virtual int GetCount();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void ZeroForce();
};

class SpindleSpecies : public Species<Spindle> {
  protected:
    bool midstep_;
  public:
    SpindleSpecies() : Species() {
      SetSID(species_id::spindle);
      midstep_ = true;
    }
    void Init(system_parameters *params, space_struct *space, long seed) {
      Species::Init(params, space, seed);
      sparams_ = &(params_->spindle);
    }
    void UpdatePositions() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        it->UpdatePosition(midstep_);
      }
      midstep_ = !midstep_;
    }
};

#endif // _SIMCORE_SPINDLE_H_

