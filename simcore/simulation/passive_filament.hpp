#ifndef _SIMCORE_PASSIVE_FILAMENT_H_
#define _SIMCORE_PASSIVE_FILAMENT_H_

#include "mesh.hpp"
#include "species.hpp"

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class PassiveFilament : public Mesh {
 private:
  int theta_validation_run_flag_, diffusion_validation_run_flag_, stoch_flag_,
      metric_forces_;
  double max_length_, min_length_, max_bond_length_, min_bond_length_,
      persistence_length_,
      friction_ratio_,  // friction_par/friction_perp
      friction_par_, friction_perp_, rand_sigma_par_, rand_sigma_perp_;

  std::vector<double> gamma_inverse_,
      tensions_,       // n_sites-1
      g_mat_lower_,    // n_sites-2
      g_mat_upper_,    // n_sites-2
      g_mat_diag_,     // n_sites-1
      det_t_mat_,      // n_sites+1
      det_b_mat_,      // n_sites+1
      g_mat_inverse_,  // n_sites-2
      k_eff_,          // n_sites-2
      h_mat_diag_,     // n_sites-1
      h_mat_upper_,    // n_sites-2
      h_mat_lower_,    // n_sites-2
      cos_thetas_;
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
  void InsertPassiveFilament(bool force_overlap = false);
  bool InsertFirstBond();
  void UpdateAvgPosition();
  void ReportAll();
  // Anchor * anchor_; //FIXME temporary?
  bool CheckBondLengths();

 public:
  PassiveFilament();
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
  double const GetPersistenceLength() { return persistence_length_; }
  double const GetBondLength() { return bond_length_; }
  int const GetNBonds() { return n_bonds_; }
  void GetAvgOrientation(double *au);
  void GetAvgPosition(double *ap);
  std::vector<double> const *const GetThetas() { return &cos_thetas_; }
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

typedef std::vector<PassiveFilament>::iterator passive_filament_iterator;
typedef std::vector<std::pair<std::vector<PassiveFilament>::iterator,
                              std::vector<PassiveFilament>::iterator> >
    passive_filament_chunk_vector;

class PassiveFilamentSpecies : public Species<PassiveFilament> {
 protected:
  bool midstep_;
  // Analysis structures
  double e_bend_, tot_angle_, mse2e_, mse2e2_;
  int **theta_histogram_;
  int time_, n_bins_, n_samples_;
  std::fstream theta_file_, mse2e_file_;

 public:
  PassiveFilamentSpecies() : Species() {
    SetSID(species_id::passive_filament);
    midstep_ = true;
  }
  void Init(system_parameters *params, space_struct *space, long seed) {
    Species::Init(params, space, seed);
    sparams_ = &(params_->passive_filament);
    if (params_->passive_filament.packing_fraction > 0) {
      if (params_->passive_filament.length <= 0) {
        error_exit(
            "Packing fraction with polydisperse lengths not implemented yet\n");
      }
      double fil_vol;
      if (params_->n_dim == 2) {
        fil_vol = params_->passive_filament.length *
                      params_->passive_filament.diameter +
                  0.25 * M_PI * SQR(params_->passive_filament.diameter);
        sparams_->num = params_->passive_filament.packing_fraction *
                        space_->volume / fil_vol;
      } else {
        fil_vol = 0.25 * M_PI * SQR(params_->passive_filament.diameter) *
                      params_->passive_filament.length +
                  M_PI * CUBE(params_->passive_filament.diameter) / 6.0;
        sparams_->num = params_->passive_filament.packing_fraction *
                        space_->volume / fil_vol;
      }
      // DPRINTF("  passive_filament_num: %d\n  sys_volume: %2.2f\n  fil_volume:
      // %2.2f\n  packing_fraction: %2.2f\n",sparams_->num, space_->volume,
      // fil_vol, params_->passive_filament.packing_fraction);
    }
  }
  void InitAnalysis();
  void InitThetaAnalysis();
  void InitMse2eAnalysis();
  void RunAnalysis();
  void RunThetaAnalysis();
  void RunMse2eAnalysis();
  void FinalizeAnalysis();
  void FinalizeMse2eAnalysis();
  void FinalizeThetaAnalysis();
  void UpdatePositions() {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    passive_filament_chunk_vector chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = members_.size() / max_threads;
    passive_filament_iterator cur_iter = members_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      passive_filament_iterator last_iter = cur_iter;
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
    for (passive_filament_iterator it = members_.begin(); it != members_.end();
         ++it)
      it->UpdatePosition(midstep_);
#endif

    midstep_ = !midstep_;
  }
  virtual void CenteredOrientedArrangement() {}
};

#endif  // _SIMCORE_PASSIVE_FILAMENT_H_
