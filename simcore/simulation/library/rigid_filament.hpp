#ifndef _SIMCORE_RIGID_FILAMENT_H_
#define _SIMCORE_RIGID_FILAMENT_H_

#include "exponential_dist.hpp"
#include "flory_schulz.hpp"
#include "mesh.hpp"

class RigidFilament : public Mesh {
 private:
  rigid_filament_parameters *sparams_;
  double gamma_par_ = 0;
  double gamma_perp_ = 0;
  double gamma_rot_ = 0;
  double diffusion_par_ = 0;
  double diffusion_perp_ = 0;
  double diffusion_rot_ = 0;
  double body_frame_[6];
  double min_length_;
  double max_length_;

  // bool force_induced_catastrophe_flag_;
  // bool diffusion_validation_run_flag_;
  bool stoch_flag_;
  // bool optical_trap_flag_;
  // bool polydispersity_flag_;

  // int optical_trap_fixed_;
  // int trapped_site_;
  int n_step_ = 0;
  int eq_steps_;
  int eq_steps_count_ = 0;

  // poly_state poly_;
  // double v_poly_;
  // double v_depoly_;
  // double p_s2g_;
  // double p_s2p_;
  // double p_p2s_;
  // double p_p2g_;
  // double p_g2s_;
  // double p_g2p_;

  // double driving_factor_;
  // double fic_factor_;
  // double tip_force_;

  // double optical_trap_spring_;
  // double optical_trap_pos_[3];
  // double optical_trap_pos2_[3];

  // double polydispersity_factor_;

  void UpdateSiteBondPositions();
  void SetDiffusion();
  void GenerateProbableOrientation();
  // void CalculateAngles();
  // void CalculateTangents();
  void AddRandomForces();

  // void ConstructUnprojectedRandomForces();  //?
  // void GeometricallyProjectRandomForces();  //?

  void UpdateSitePositions();

  void ApplyForcesTorques();
  void ApplyInteractionForces();

  void SetParameters();
  void InitRigidFilamentLength();
  void UpdateAvgPosition();
  void DynamicInstability();
  void UpdatePolyState();
  void GrowFilament();
  void ReportAll();

  void CalculateBinding();
  // bool CheckBondLengths();

 protected:
  void InsertRigidFilament(std::string insertion_type, double buffer = -1);
  void GetBodyFrame();
  void AddRandomDisplacement();
  void AddRandomReorientation();

 public:
  RigidFilament(unsigned long seed);
  virtual void Init(rigid_filament_parameters *sparams);
  virtual void InsertAt(const double *const new_pos, const double *const u);
  virtual void Integrate();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void UpdatePosition();
  double const GetLength() { return length_; }
  // double const GetDriving() { return driving_factor_; }
  // void Reserve();
  void CheckFlocking();

  void GetNematicOrder(double *nematic_order_tensor);
  void GetPolarOrder(double *polar_order_vector);
  double GetTipZ() { return sites_[n_sites_ - 1].GetOrientation()[n_dim_ - 1]; }
  double const *const GetHeadPosition() {
    return sites_[n_sites_ - 1].GetPosition();
  }
  // double const *const GetTailPosition() { return sites_[0].GetPosition(); }
  // double const *const GetTailOrientation() {
  //  return sites_[0].GetOrientation();
  //}
  void WritePosit(std::fstream &oposit);
  void ReadPosit(std::fstream &iposit);
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadCheckpoint(std::fstream &icheck);
  void ScalePosition();
  double const GetVolume();
};

typedef std::vector<RigidFilament>::iterator rigid_filament_iterator;
typedef std::vector<std::pair<std::vector<RigidFilament>::iterator,
                              std::vector<RigidFilament>::iterator>>
    rigid_filament_chunk_vector;

#endif  // _SIMCORE_RIGID_FILAMENT_H_
