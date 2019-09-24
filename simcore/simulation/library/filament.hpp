#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "mesh.hpp"
#include "species.hpp"
//#include "motor.hpp"

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif
#include "flory_schulz.hpp"
#include "minimum_distance.hpp"

class Filament : public Mesh {
private:
  int dynamic_instability_flag_;
  int force_induced_catastrophe_flag_;
  int theta_validation_run_flag_;
  int diffusion_validation_run_flag_;
  int spiral_flag_;
  int shuffle_flag_;
  int stoch_flag_;
  int flagella_flag_;
  int metric_forces_;
  int optical_trap_flag_;
  int cilia_trap_flag_;
  int optical_trap_fixed_;
  int trapped_site_;
  // TEMPORARY FIXME
  int n_step_ = 0;
  int n_motors_bound_;
  int eq_steps_;
  int eq_steps_count_ = 0;
  int polydispersity_flag_ = 0;
  int polydispersity_warn_on_truncate_ = 0;
  double min_length_;
  double max_bond_length_;
  double min_bond_length_;
  double persistence_length_;
  double persistence_length_target_;
  double target_step_;
  double friction_ratio_; // friction_par/friction_perp
  double friction_par_;
  double friction_perp_;
  double rand_sigma_par_;
  double rand_sigma_perp_;
  double flagella_freq_;
  double flagella_period_;
  double flagella_amplitude_;
  double v_poly_;
  double v_depoly_;
  double p_s2g_;
  double p_s2p_;
  double p_p2s_;
  double p_p2g_;
  double p_g2s_;
  double p_g2p_;
  double driving_factor_;
  double fic_factor_;
  double curvature_ = 0;
  double shuffle_factor_;
  double shuffle_frequency_;
  double spiral_number_;
  double tip_force_;
  double optical_trap_spring_;
  double optical_trap_pos_[3];
  double optical_trap_pos2_[3];
  double max_length_;
  double polydispersity_factor_;
  std::vector<double> gamma_inverse_;
  std::vector<double> tensions_;      // n_sites-1
  std::vector<double> g_mat_lower_;   // n_sites-2
  std::vector<double> g_mat_upper_;   // n_sites-2
  std::vector<double> g_mat_diag_;    // n_sites-1
  std::vector<double> det_t_mat_;     // n_sites+1
  std::vector<double> det_b_mat_;     // n_sites+1
  std::vector<double> g_mat_inverse_; // n_sites-2
  std::vector<double> k_eff_;         // n_sites-2
  std::vector<double> h_mat_diag_;    // n_sites-1
  std::vector<double> h_mat_upper_;   // n_sites-2
  std::vector<double> h_mat_lower_;   // n_sites-2
  std::vector<double> cos_thetas_;
  poly_state poly_;
  void UpdateSiteBondPositions();
  void SetDiffusion();
  void GenerateProbableOrientation();
  void CalculateAngles(bool rescale = true);
  void CalculateTangents();
  void AddRandomForces();
  void ConstructUnprojectedRandomForces();
  void GeometricallyProjectRandomForces();
  void CalculateBendingForces();
  void CalculateTensions();
  void UpdateSitePositions();
  void UpdateSiteOrientations();
  void ApplyForcesTorques();
  // void ApplyAnchorForces(); // FIXME temporary
  void ApplyInteractionForces();
  void SetParameters();
  void InitElements();
  void InsertFilament(bool force_overlap = false);
  bool InsertFirstBond();
  void UpdateAvgPosition();
  void DynamicInstability();
  void UpdatePolyState();
  void GrowFilament();
  void RescaleBonds();
  void InitSpiral2D();
  void ReportAll();
  // std::vector<Motor> motors_; //FIXME temporary
  // Anchor * anchor_; //FIXME temporary?
  // void UnbindMotor();
  // void BindMotor();
  void CalculateBinding();
  // void RebindMotors();
  bool CheckBondLengths();

public:
  Filament();
  virtual void Init(bool force_overlap = false);
  virtual void InsertAt(double *pos, double *u);
  // virtual void SetAnchor(Anchor * a);
  virtual bool CheckBounds(double buffer = 0);
  // void DiffusionValidationInit();
  virtual void Integrate();
  virtual void Draw(std::vector<graph_struct *> *graph_array);
  virtual void UpdatePosition() {}
  virtual void UpdatePosition(bool midstep);
  double const GetLength() { return length_; }
  double const GetDriving() { return driving_factor_; }
  double const GetPersistenceLength() { return persistence_length_; }
  void CheckFlocking();
  int const GetNBonds() { return n_bonds_; }
  std::vector<double> const *const GetThetas() { return &cos_thetas_; }
  void CalculateSpiralNumber();
  double GetSpiralNumber();
  void GetNematicOrder(double *nematic_order_tensor);
  void GetPolarOrder(double *polar_order_vector);
  double GetTipZ() { return sites_[n_sites_ - 1].GetOrientation()[n_dim_ - 1]; }
  double const *const GetHeadPosition() {
    return sites_[n_sites_ - 1].GetPosition();
  }
  double const *const GetTailPosition() { return sites_[0].GetPosition(); }
  double const *const GetTailOrientation() {
    return sites_[0].GetOrientation();
  }
  void AddTorqueTail(double *t) { bonds_[0].AddTorque(t); }
  void AddForceTail(double *f) { sites_[0].AddForce(f); }
  void WritePosit(std::fstream &oposit);
  void ReadPosit(std::fstream &iposit);
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadCheckpoint(std::fstream &icheck);
  void ScalePosition();
  double const GetVolume();
};

typedef std::vector<Filament>::iterator filament_iterator;
typedef std::vector<
    std::pair<std::vector<Filament>::iterator, std::vector<Filament>::iterator>>
    filament_chunk_vector;

class FilamentSpecies : public Species<Filament> {
protected:
  bool midstep_;
  // Analysis structures
  MinimumDistance mindist_;
  double e_bend_, tot_angle_, mse2e_, mse2e2_, polar_bin_width_,
      contact_bin_width_, contact_cut_, po_avg_, cn_avg_,
      *nematic_order_tensor_, *polar_order_vector_;
  int **theta_histogram_;
  int *polar_order_histogram_;
  int time_, n_bins_, n_bins_1d_, n_samples_, orientation_corr_n_steps_;
  std::fstream spiral_file_, flock_file_, theta_file_, mse2e_file_,
      global_order_file_, local_order_file_, polar_order_file_,
      orientation_corr_file_, crossing_file_, polar_order_avg_file_,
      in_out_file_;

public:
  FilamentSpecies() : Species() {
    SetSID(species_id::filament);
    midstep_ = true;
  }
  void Init(system_parameters *params, space_struct *space, long seed) {
    Species::Init(params, space, seed);
    sparams_ = &(params_->filament);
    if (params_->filament.packing_fraction > 0) {
      if (params_->filament.length <= 0) {
        error_exit(
            "Packing fraction with polydisperse lengths not implemented yet\n");
      }
      double fil_vol;
      if (params_->n_dim == 2) {
        fil_vol = params_->filament.length * params_->filament.diameter +
                  0.25 * M_PI * SQR(params_->filament.diameter);
        sparams_->num =
            params_->filament.packing_fraction * space_->volume / fil_vol;
      } else {
        fil_vol = 0.25 * M_PI * SQR(params_->filament.diameter) *
                      params_->filament.length +
                  M_PI * CUBE(params_->filament.diameter) / 6.0;
        sparams_->num =
            params_->filament.packing_fraction * space_->volume / fil_vol;
      }
      // DPRINTF("  filament_num: %d\n  sys_volume: %2.2f\n  fil_volume: %2.2f\n
      // packing_fraction: %2.2f\n",sparams_->num, space_->volume, fil_vol,
      // params_->filament.packing_fraction);
    }
  }

  void InitAnalysis();
  void RunAnalysis();
  void FinalizeAnalysis();

  void InitSpiralAnalysis();
  void RunSpiralAnalysis();
  // void FinalizeSpiralAnalysis();

  void InitCrossingAnalysis();
  void RunCrossingAnalysis();
  void FinalizeCrossingAnalysis();

  void InitThetaAnalysis();
  void RunThetaAnalysis();
  void FinalizeThetaAnalysis();

  void InitMse2eAnalysis();
  void RunMse2eAnalysis();
  void FinalizeMse2eAnalysis();

  void InitGlobalOrderAnalysis();
  void RunGlobalOrderAnalysis();
  void FinalizeGlobalOrderAnalysis();

  void InitPolarOrderAnalysis();
  void RunPolarOrderAnalysis();
  void FinalizePolarOrderAnalysis();

  void InitLocalOrderAnalysis();
  void RunLocalOrderAnalysis();
  void WriteLocalOrderData();
  void FinalizeLocalOrderAnalysis();

  void InitOrientationCorrelationAnalysis();
  void RunOrientationCorrelationAnalysis();
  void FinalizeOrientationCorrelationAnalysis();

  void InitFlockingAnalysis();
  void RunFlockingAnalysis();
  void FinalizeFlockingAnalysis();

  void UpdatePositions() {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    filament_chunk_vector chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = members_.size() / max_threads;
    filament_iterator cur_iter = members_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      filament_iterator last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
    {
#pragma omp for
      for (int i = 0; i < max_threads; ++i)
        for (auto it = chunks[i].first; it != chunks[i].second; ++it)
          it->UpdatePosition(midstep_);
    }
#else
    for (filament_iterator it = members_.begin(); it != members_.end(); ++it)
      it->UpdatePosition(midstep_);
#endif

    midstep_ = !midstep_;
  }
  virtual void CenteredOrientedArrangement() {}
};

#endif // _SIMCORE_FILAMENT_H_
