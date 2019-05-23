#ifndef _SIMCORE_SPHEROCYLINDER_H_
#define _SIMCORE_SPHEROCYLINDER_H_

#include "species.h"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class Spherocylinder : public Object {
  protected:
    double gamma_par_,
           gamma_perp_,
           gamma_rot_,
           diffusion_par_,
           diffusion_perp_,
           diffusion_rot_,
           body_frame_[6];
    bool is_midstep_;
    void ApplyForcesTorques();
    void InsertSpherocylinder();
    void SetDiffusion();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void AddRandomReorientation();
    void Integrate();
  public:
    Spherocylinder();
    void Init();
    void UpdatePosition();
};

class SpherocylinderSpecies : public Species<Spherocylinder> {
  protected:
    bool midstep_;
    double ** pos0_,
           ** u0_,
           * msd_,
           * msd_err_,
           * vcf_,
           * vcf_err_;
    int time_,
        time_avg_interval_,
        n_samples_;
    std::fstream diff_file_;
    void InitDiffusionAnalysis();
    void DiffusionAnalysis();
    void CalculateMSD();
    void CalculateVCF();
    void UpdateInitPositions();
    void FinalizeDiffusionAnalysis();
  public:
    SpherocylinderSpecies() : Species() {
      SetSID(species_id::spherocylinder);
    }
    void Init(system_parameters *params, space_struct *space, long seed) {
      Species::Init(params, space, seed);
      sparams_ = &(params_->spherocylinder);
      midstep_ = params_->spherocylinder.midstep;
    }
    void UpdatePositions() {
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        it->UpdatePosition();
      }
    }
    virtual void InitAnalysis();
    virtual void RunAnalysis();
    virtual void FinalizeAnalysis();
};

#endif // _SIMCORE_SPHEROCYLINDER_H_

