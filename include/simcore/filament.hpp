#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "exponential_dist.hpp"
#include "flory_schulz.hpp"
#include "mesh.hpp"

class Filament : public Mesh {
 private:
  filament_parameters *sparams_;
  bool force_induced_catastrophe_flag_ = false;
  bool theta_validation_run_flag_ = false;
  bool diffusion_validation_run_flag_ = false;
  bool spiral_flag_ = false;
  bool stoch_flag_ = true;
  bool flagella_flag_ = false;
  bool optical_trap_flag_ = false;
  bool cilia_trap_flag_ = false;
  bool normalize_switch_ = true;
  bool nematic_driving_ = false;
  bool custom_set_tail_ = false;
  int n_normalize_ = 0;
  int optical_trap_fixed_ = 0;
  int trapped_site_ = 0;
  int n_step_ = 0;
  int eq_steps_ = 0;
  int in_flock_ = 0;           // 0 if not in flock, 1 if interior, 2 if exterior
  int flock_change_state_ = 0; // 0 if same as previous step, 1 if joined flock, 2
                               // if left flock
  double min_length_;
  double max_length_;
  double max_bond_length_;
  double min_bond_length_;
  double bending_stiffness_;
  double friction_ratio_;  // friction_par/friction_perp
  double friction_par_;
  double friction_perp_;
  double rand_sigma_par_;
  double rand_sigma_perp_;
  double flagella_freq_;
  double flagella_period_;
  double flagella_amplitude_;
  double p_driving_switch_ = 0;
  double v_poly_;
  double v_depoly_;
  double p_s2g_;
  double p_s2p_;
  double p_p2s_;
  double p_p2g_;
  double p_g2s_;
  double p_g2p_;
  double driving_factor_ = 0;
  double peclet_number_;
  double flexure_number_;
  double fic_factor_;
  double curvature_ = 0;
  double spiral_number_;
  double tip_force_;
  double optical_trap_spring_;
  double optical_trap_pos_[3];
  double optical_trap_pos2_[3];
  double polydispersity_factor_;
  bool error_analysis_;
  std::vector<int> error_rates_;
  std::vector<double> gamma_inverse_;
  std::vector<double> tensions_;       // n_sites-1
  std::vector<double> g_mat_lower_;    // n_sites-2
  std::vector<double> g_mat_upper_;    // n_sites-2
  std::vector<double> g_mat_diag_;     // n_sites-1
  std::vector<double> det_t_mat_;      // n_sites+1
  std::vector<double> det_b_mat_;      // n_sites+1
  std::vector<double> g_mat_inverse_;  // n_sites-2
  std::vector<double> k_eff_;          // n_sites-2
  std::vector<double> h_mat_diag_;     // n_sites-1
  std::vector<double> h_mat_upper_;    // n_sites-2
  std::vector<double> h_mat_lower_;    // n_sites-2
  std::vector<double> cos_thetas_;
  poly_state poly_;
  void UpdateSiteBondPositions();
  void SetDiffusion();
  void GenerateProbableOrientation();
  void CalculateAngles();
  void CalculateTangents();
  void AddRandomForces();
  void ConstructUnprojectedRandomForces();
  void GeometricallyProjectRandomForces();
  void CalculateBendingForces();
  void CalculateTensions();
  void UpdateSitePositions();
  void RotateToReferenceFrame();
  void ApplyForcesTorques();
  void ApplyInteractionForces();
  void SetParameters();
  void InitFilamentLength();
  void InsertFilament();
  void InsertFirstBond();
  void UpdateAvgPosition();
  void DynamicInstability();
  void UpdatePolyState();
  void GrowFilament();
  void RescaleBonds();
  void InitSpiral2D();
  void ReportAll();
  void CalculateBinding();
  bool CheckBondLengths();
 public:
  Filament(unsigned long seed);
  virtual void Init(filament_parameters *sparams);
  virtual void InsertAt(const double *const new_pos, const double *const u);
  virtual void Integrate();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void UpdatePosition() {}
  virtual void UpdatePosition(bool midstep);
  double const GetLength() { return length_; }
  double const GetDriving() { return driving_factor_; }
  double const GetPersistenceLength() { return bending_stiffness_; }
  void Reserve();
  void CheckFlocking();
  const int GetNBonds() { return n_bonds_; }
  const int GetHandedness() { return SIGNOF(curvature_); }
  const int GetFlockChangeState() { return flock_change_state_; }
  const int GetFlockType() { return in_flock_; }
  std::vector<double> const *const GetThetas() { return &cos_thetas_; }
  void CalculateSpiralNumber();
  double GetSpiralNumber();
  const double GetCenterOfCurvature(double *center);
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
  void GetErrorRates(std::vector<int> &rates);
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
  // const double GetLength() { return length_; };
};

typedef std::vector<Filament>::iterator filament_iterator;
typedef std::vector<
    std::pair<std::vector<Filament>::iterator, std::vector<Filament>::iterator>>
    filament_chunk_vector;

#endif  // _SIMCORE_FILAMENT_H_
