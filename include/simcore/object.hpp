#ifndef _SIMCORE_OBJECT_H_
#define _SIMCORE_OBJECT_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "rng.hpp"
#include <mutex>

class Object {
private:
  int oid_;
  int mesh_id_;
  static int _next_oid_;
  static std::mutex _obj_mtx_;
  void InitOID();

protected:
  static system_parameters *params_;
  static space_struct *space_;
  static int n_dim_;
  static double delta_;
  species_id sid_;
  obj_type type_;
  graph_struct g_;
  RNG rng_;
  draw_type draw_;
  int n_contact_;
  int in_flock_;           // 0 if not in flock, 1 if interior, 2 if exterior
  int flock_change_state_; // 0 if same as previous step, 1 if joined flock, 2
                           // if left flock
  double position_[3];
  double prev_position_[3];
  double prev_orientation_[3];
  double scaled_position_[3];
  double orientation_[3];
  double force_[3];
  double torque_[3];
  double dr_zero_[3];
  double color_;
  double diameter_;
  double length_;
  double p_energy_;
  double dr_tot_;
  double polar_order_;
  double contact_number_;
  bool interacting_;
  bool is_mesh_;
  bool has_overlap_;
  bool interactor_update_;

  std::vector<Object *> interactors_;
  std::vector<object_interaction> ixs_;
  void UpdateKMC();

public:
  Object(unsigned long seed);
  virtual ~Object() = default;
  // kmcx parameter
  int gid;
  double length;
  double radius;
  double pos[3];
  double direction[3];

  // Static functions
  static void SetParams(system_parameters *params);
  static void SetSpace(space_struct *space);
  static void SetNDim(int n_dim);
  static void SetDelta(double delta);
  static const double GetDelta();
  static const int GetNextOID();
  static void SetNextOID(const int next_oid);
  // Trivial Set/Get functions
  void SetSID(species_id sid);
  void SetType(obj_type type);
  void SetPosition(const double *const new_pos);
  void SetScaledPosition(const double *const spos);
  void SetOrientation(const double *const u);
  void SetPrevPosition(const double *const ppos);
  void SetPrevOrientation(const double *const pu);
  void SetDiameter(double new_diameter);
  void SetLength(double new_length);
  void AddForce(const double *const f);
  void SubForce(const double *const f);
  void SetForce(const double *const f);
  void AddTorque(const double *const t);
  void SubTorque(const double *const t);
  void SetTorque(const double *const t);
  void AddPotential(const double p);
  void AddPolarOrder(const double po);
  void AddContactNumber(const double cn);
  void SetInteractor(bool ix);
  void ToggleIsMesh();
  void CalcPolarOrder();
  void ZeroPolarOrder();
  species_id const GetSID();
  obj_type const GetType();
  const int GetOID() const;
  const int GetMeshID() const;
  const double *const GetPosition();
  const double *const GetPrevPosition();
  const double *const GetPrevOrientation();
  const double *const GetScaledPosition();
  const double *const GetOrientation();
  virtual void GetAvgPosition(double *ap);
  virtual void GetAvgOrientation(double *au);
  virtual void SetAvgPosition();
  const double GetDiameter();
  const double GetLength();
  const double *const GetForce();
  const double *const GetTorque();
  const double GetPotentialEnergy();
  const double GetPolarOrder();
  const double GetContactNumber();
  const bool IsInteractor();
  const bool IsMesh();
  const bool CheckInteractorUpdate();
  void HasOverlap(bool overlap);
  void SetFlockType(int in_flock);
  void SetFlockChangeState(int fcs);
  int GetFlockType();
  int GetFlockChangeState();
  void SetOID(int oid);
  void SetMeshID(int mid);

  // Virtual functions
  virtual void Init(species_base_parameters *sparams) {}
  virtual void InsertRandom(double buffer = -1);
  virtual void InsertRandomOriented(const double *const u);
  virtual void InsertAt(const double *const new_pos, const double *const u);
  virtual void ZeroForce();
  virtual void UpdatePeriodic();
  virtual void UpdatePosition() {}
  virtual void ResetPreviousPosition();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void SetColor(const double c, draw_type dtype);
  virtual void ScalePosition();
  virtual int GetCount();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual const double *const GetInteractorPosition();
  virtual const double *const GetInteractorPrevPosition();
  virtual const double *const GetInteractorScaledPosition();
  virtual const double *const GetInteractorOrientation();
  virtual const double GetInteractorDiameter();
  virtual const double GetInteractorLength();
  virtual const double GetVolume();
  virtual void UpdateDrTot();
  virtual const double GetDrTot();
  virtual void ZeroDrTot();
  virtual bool HasNeighbor(int other_id);
  virtual void GiveInteraction(object_interaction ix);
  virtual void ApplyInteractions();
  virtual void FlagDuplicateInteractions();
  // virtual std::vector<Interaction *> *GetInteractions();
  virtual void ClearInteractions();
  virtual void Cleanup();
  // virtual void BindAnchor(anchor *ix);
  // virtual void UnbindAnchor();

  // I/O functions
  virtual void Report();
  virtual void WritePosit(std::fstream &oposit);
  virtual void ReadPosit(std::fstream &iposit);
  virtual void WriteSpec(std::fstream &ospec);
  virtual void ReadSpec(std::fstream &ispec);
  virtual void ReadPositFromSpec(std::fstream &ispec);
  virtual void WriteCheckpoint(std::fstream &ocheck);
  virtual void WriteCheckpointHeader(std::fstream &ocheck);
  virtual void ReadCheckpoint(std::fstream &icheck);
  virtual void ReadCheckpointHeader(std::fstream &icheck);
};

// void MinimumDistance(Object* o1, Object* o2, Interaction *ix, space_struct
// *space); void BoundaryConditions(Object * o1, space_struct *space);

#endif // _SIMCORE_OBJECT_H_
