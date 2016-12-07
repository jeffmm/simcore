#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

class Filament : public Composite<Site,Bond> {

  private:
    int n_bonds_,
        n_sites_,
        dynamic_instability_flag_,
        force_induced_catastrophe_flag_,
        theta_validation_flag_,
        diffusion_validation_flag_,
        metric_forces_;
    double max_length_,
           min_length_,
           max_child_length_,
           child_length_,
           persistence_length_,
           gamma_ratio_, // gamma_par/gamma_perp
           gamma_par_,
           gamma_perp_,
           rand_sigma_par_,
           rand_sigma_perp_,
           v_poly_,
           v_depoly_,
           p_s2g_,
           p_s2p_,
           p_p2s_,
           p_p2g_,
           p_g2s_,
           p_g2p_,
           driving_factor_,
           tip_force_;
    std::vector<double> gamma_inverse_,
                        tensions_, //n_sites-1
                        g_mat_lower_, //n_sites-2
                        g_mat_upper_, //n_sites-2
                        g_mat_diag_, //n_sites-1
                        det_t_mat_, //n_sites+1
                        det_b_mat_, //n_sites+1
                        g_mat_inverse_, //n_sites-2
                        k_eff_, //n_sites-2
                        h_mat_diag_, //n_sites-1
                        h_mat_upper_, //n_sites-2
                        h_mat_lower_, //n_sites-2
                        cos_thetas_;
    poly_state_t poly_state_;
    void UpdateSiteBondPositions();
    void SetDiffusion();
    void GenerateProbableOrientation();
    void CalculateAngles();
    void CalculateTangents();
    void UpdatePrevPositions();
    void AddRandomForces();
    void ConstructUnprojectedRandomForces();
    void GeometricallyProjectRandomForces();
    void CalculateBendingForces();
    void CalculateTensions();
    void UpdateSitePositions(bool midstep);
    void UpdateBondPositions();
    void ApplyForcesTorques();
    void SetParameters(system_parameters *params);
    void InitElements(system_parameters *params, space_struct *space);
    void UpdateAvgPosition();
    //void DynamicInstability();
    void DumpAll();

  public:
    Filament(system_parameters *params, space_struct * space, 
        long seed, SID sid) : Composite(params, space, seed, sid) {
      SetParameters(params);
      InitElements(params, space);
    }
    virtual void Init();
    void DiffusionValidationInit();
    virtual void Integrate(bool midstep);
    virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition() {}
    virtual void UpdatePosition(bool midstep);
    void GetAvgOrientation(double * au);
    void GetAvgPosition(double * ap);
    std::vector<double> const * const GetThetas() {
      return &cos_thetas_;
    }
    void WritePosit(std::fstream &op);
};

class FilamentSpecies : public Species<Filament> {
  protected:
    //void InitPotentials(system_parameters *params);
    //DiffusionProperties diffusion_properties_;
    bool theta_validation_,
         diffusion_validation_,
         midstep_;
    int ***theta_distribution_;
    int nbins_,
        ibin_,
        n_dim_,
        nvalidate_,
        ivalidate_;
    void ValidateDiffusion();
    void ValidateThetaDistributions();
    void WriteThetaValidation(std::string run_name);
    void WriteDiffusionValidation(std::string run_name);
  public:
    FilamentSpecies() : Species() {
      SetSID(SID::filament);
      midstep_ = true;
      ibin_ = 0;
    }
    ~FilamentSpecies() {
      if (theta_validation_) {
        for (int i=0; i<n_members_; ++i) {
          for (int j=0; j<7; ++j) {
            delete[] theta_distribution_[i][j];
          }
          delete[] theta_distribution_[i];
        }
        delete[] theta_distribution_;
      }
    }
    FilamentSpecies(const FilamentSpecies& that) {
      theta_distribution_ = that.theta_distribution_;
      theta_validation_ = that.theta_validation_;
      midstep_ = that.midstep_;
      ibin_ = that.ibin_;
      nbins_ = that.nbins_;
    }
    FilamentSpecies& operator=(FilamentSpecies const& that) {
      theta_distribution_ = that.theta_distribution_;
      theta_validation_ = that.theta_validation_;
      midstep_ = that.midstep_;
      ibin_ = that.ibin_;
      nbins_ = that.nbins_;
      return *this;
    }
    void Init() {
      Species::Init();
    }
    void UpdatePositions() {
      if (diffusion_validation_ && ivalidate_%nvalidate_ == 0)
        ValidateDiffusion();
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->UpdatePosition(midstep_);
      }
      if (theta_validation_ && midstep_)
        ValidateThetaDistributions();
      midstep_ = !midstep_;
      ivalidate_++;
    }
    void WriteOutputs(std::string run_name);
    void Configurator();

};

#endif // _SIMCORE_FILAMENT_H_
