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
    void GenerateRandomForces();
    void ProjectRandomForces();
    void CalculateBendingForces();
    void CalculateTensions();
    void UpdateSitePositions(bool midstep);
    void UpdateBondPositions();
    void ApplyForcesTorques();
    void SetParameters(system_parameters *params);
    void InitElements(system_parameters *params, space_struct *space);
    //void UpdateOrientation();
    //void AddRandomDisplacement();
    //void ApplyForcesTorques();
    //void DynamicInstability();
    void DumpAll();

  public:
    Filament(system_parameters *params, space_struct * space, 
        long seed, SID sid) : Composite(params, space, seed, sid) {
      SetParameters(params);
      InitElements(params, space);
    }
    ~Filament() {}
    Filament(const Filament& that) : Composite(that) {
      n_bonds_=that.n_bonds_;
      n_sites_ = that.n_sites_;
      dynamic_instability_flag_ = that.dynamic_instability_flag_;
      force_induced_catastrophe_flag_ = that.force_induced_catastrophe_flag_;
      metric_forces_ = that.metric_forces_;
      theta_validation_flag_ = that.theta_validation_flag_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_ratio_ = that.gamma_ratio_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      v_poly_ = that.v_poly_;
      v_depoly_ = that.v_depoly_;
      p_s2g_ = that.p_s2g_;
      p_s2p_ = that.p_s2p_;
      p_p2s_ = that.p_p2s_;
      p_p2g_ = that.p_p2g_;
      p_g2s_ = that.p_g2s_;
      p_g2p_ = that.p_g2p_;
      tip_force_ = that.tip_force_;
      gamma_inverse_ = that.gamma_inverse_;
      tensions_ = that.tensions_;
      g_mat_lower_ = that.g_mat_lower_;
      g_mat_upper_ = that.g_mat_upper_;
      g_mat_diag_ = that.g_mat_diag_;
      h_mat_diag_ = that.h_mat_diag_;
      h_mat_lower_ = that.h_mat_lower_;
      h_mat_upper_ = that.h_mat_upper_;
      cos_thetas_ = that.cos_thetas_;
      k_eff_ = that.k_eff_;
      g_mat_inverse_ = that.g_mat_inverse_;
      det_t_mat_ = that.det_t_mat_;
      det_b_mat_ = that.det_b_mat_;
      poly_state_ = that.poly_state_;
    }
    Filament& operator=(Filament const& that) {
      Composite::operator=(that); 
      n_bonds_=that.n_bonds_;
      n_sites_ = that.n_sites_;
      dynamic_instability_flag_ = that.dynamic_instability_flag_;
      force_induced_catastrophe_flag_ = that.force_induced_catastrophe_flag_;
      theta_validation_flag_ = that.theta_validation_flag_;
      metric_forces_ = that.metric_forces_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_ratio_ = that.gamma_ratio_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      v_poly_ = that.v_poly_;
      v_depoly_ = that.v_depoly_;
      p_s2g_ = that.p_s2g_;
      p_s2p_ = that.p_s2p_;
      p_p2s_ = that.p_p2s_;
      p_p2g_ = that.p_p2g_;
      p_g2s_ = that.p_g2s_;
      p_g2p_ = that.p_g2p_;
      tip_force_ = that.tip_force_;
      gamma_inverse_ = that.gamma_inverse_;
      tensions_ = that.tensions_;
      g_mat_lower_ = that.g_mat_lower_;
      g_mat_upper_ = that.g_mat_upper_;
      g_mat_diag_ = that.g_mat_diag_;
      h_mat_diag_ = that.h_mat_diag_;
      h_mat_lower_ = that.h_mat_lower_;
      h_mat_upper_ = that.h_mat_upper_;
      cos_thetas_ = that.cos_thetas_;
      k_eff_ = that.k_eff_;
      g_mat_inverse_ = that.g_mat_inverse_;
      det_t_mat_ = that.det_t_mat_;
      det_b_mat_ = that.det_b_mat_;
      poly_state_ = that.poly_state_;
      return *this;
    } 
    virtual void Init();
    void DiffusionInit();
    virtual void Integrate(bool midstep);
    virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition() {}
    virtual void UpdatePositionMP() {}
    virtual void UpdatePosition(bool midstep);
    virtual void UpdatePositionMP(bool midstep);
    void GetAvgOrientation(double * au);
    void GetAvgPosition(double * ap);
    std::vector<double> const * const GetThetas() {
      return &cos_thetas_;
    }
};

class FilamentSpecies : public Species<Filament> {
  protected:
    //void InitPotentials(system_parameters *params);
    //DiffusionProperties diffusion_properties_;
    bool theta_validation_,
         diffusion_validation_,
         midstep_;
    int ***theta_distribution_;
    double ***orientations_;
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
      if (diffusion_validation_) {
        for (int i=0; i<n_members_; ++i) {
          for (int j=0; j<2; ++j) {
            delete[] orientations_[i][j];
          }
          delete[] orientations_[i];
        }
        delete[] orientations_;
      }
    }
    FilamentSpecies(const FilamentSpecies& that) : Species(that) {
      theta_distribution_ = that.theta_distribution_;
      theta_validation_ = that.theta_validation_;
      midstep_ = that.midstep_;
      ibin_ = that.ibin_;
      nbins_ = that.nbins_;
    }
    FilamentSpecies& operator=(FilamentSpecies const& that) {
      Species::operator=(that);
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
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->UpdatePosition(midstep_);
      }
    }
    void UpdatePositionsMP() {
      if (diffusion_validation_ && ivalidate_%nvalidate_ == 0)
        ValidateDiffusion();
      for (auto it=members_.begin(); it!=members_.end(); ++it) {
        (*it)->UpdatePositionMP(midstep_);
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
