#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "species.h"
#include "mesh.h"

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class Filament : public Mesh {

  private:
    int dynamic_instability_flag_,
        force_induced_catastrophe_flag_,
        theta_validation_flag_,
        diffusion_validation_flag_,
        spiral_flag_,
        stoch_flag_,
        metric_forces_,
        eq_steps_,
        eq_steps_count_ = 0;
    double max_length_,
           min_length_,
           max_bond_length_,
           persistence_length_,
           friction_ratio_, // friction_par/friction_perp
           friction_par_,
           friction_perp_,
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
    //poly_state poly_;
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
    void UpdateSiteOrientations();
    void ApplyForcesTorques();
    void SetParameters();
    void InitElements();
    void InitRandom();
    void UpdateAvgPosition();
    //void InitSpiral2D();
    void ReportAll();

  public:
    Filament();
    virtual void Init();
    virtual void InsertAt(double *pos, double *u);
    //void DiffusionValidationInit();
    virtual void Integrate(bool midstep);
    virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition() {}
    virtual void UpdatePosition(bool midstep);
    double const GetLength() { return length_;}
    double const GetDriving() {return driving_factor_;}
    double const GetPersistenceLength() {return persistence_length_;}
    double const GetBondLength() {return bond_length_;}
    int const GetNBonds() {return n_bonds_;}
    void GetAvgOrientation(double * au);
    void GetAvgPosition(double * ap);
    std::vector<double> const * const GetThetas() {
      return &cos_thetas_;
    }
    double GetTipZ() {
      return sites_[n_sites_-1].GetOrientation()[n_dim_-1];
    }
    double const * const GetHeadPos() {
      return sites_[n_sites_-1].GetPosition();
    }
    double const * const GetTailPos() {
      return sites_[0].GetPosition();
    }
    void WritePosit(std::fstream &oposit);
    void ReadPosit(std::fstream &iposit);
    void WriteSpec(std::fstream &ospec);
    void ReadSpec(std::fstream &ispec);
    void WriteCheckpoint(std::fstream &ocheck);
    void ReadCheckpoint(std::fstream &icheck);
    void ScalePosition();
};

class FilamentSpecies : public Species<Filament> {
  protected:
    bool midstep_;
    // Analysis structures
    double e_bend_,
           tot_angle_,
           mse2e_,
           mse2e2_;
    int **theta_histogram_;
    int time_,
        n_bins_,
        n_samples_;
    std::fstream spiral_file_,
                 theta_file_,
                 mse2e_file_;
  public:
    FilamentSpecies() : Species() {
      SetSID(species_id::filament);
      midstep_ = true;
    }
    void Init(system_parameters *params, space_struct *space, long seed) {
      Species::Init(params, space, seed);
      sparams_ = &(params_->filament);
    }
    void InitAnalysis();
    void InitSpiralAnalysis();
    void InitThetaAnalysis();
    void InitMse2eAnalysis();
    void RunAnalysis();
    void RunSpiralAnalysis();
    void RunThetaAnalysis();
    void RunMse2eAnalysis();
    void FinalizeAnalysis();
    void FinalizeMse2eAnalysis();
    void FinalizeThetaAnalysis();
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
            it->UpdatePosition(midstep_);
      }
#else
      for (auto it=members_.begin(); it!=members_.end(); ++it) 
        it->UpdatePosition(midstep_);
#endif

      midstep_ = !midstep_;
    }
};

#endif // _SIMCORE_FILAMENT_H_
