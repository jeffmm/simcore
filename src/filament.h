#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include <yaml-cpp/yaml.h>
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

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
    void ReadPosit(std::fstream &ip);
    void ScalePosition();
};

class FilamentSpecies : public Species<Filament> {
  protected:
    bool theta_validation_,
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
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      std::vector<std::pair<std::vector<Filament*>::iterator, std::vector<Filament*>::iterator> > chunks;
      chunks.reserve(max_threads); 
      size_t chunk_size= members_.size() / max_threads;
      auto cur_iter = members_.begin();
      for(int i = 0; i < max_threads - 1; ++i) {
        auto last_iter = cur_iter;
        std::advance(cur_iter, chunk_size);
        chunks.push_back(std::make_pair(last_iter, cur_iter));
      }
      chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
      {
#pragma omp for 
        for(int i = 0; i < max_threads; ++i)
          for(auto it = chunks[i].first; it != chunks[i].second; ++it)
            (*it)->UpdatePosition(midstep_);
      }
#else
      for (auto it=members_.begin(); it!=members_.end(); ++it) 
        (*it)->UpdatePosition(midstep_);
#endif

      //if (theta_validation_ && midstep_)
        //ValidateThetaDistributions();
      midstep_ = !midstep_;
      //ivalidate_++;
    }
    void WriteOutputs(std::string run_name);
    void Configurator();

};

#endif // _SIMCORE_FILAMENT_H_
